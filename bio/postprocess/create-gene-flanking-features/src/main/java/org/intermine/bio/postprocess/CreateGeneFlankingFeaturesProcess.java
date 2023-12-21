package org.intermine.bio.postprocess;

/*
 * Copyright (C) 2002-2022 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import org.apache.log4j.Logger;

import org.intermine.bio.util.BioQueries;
import org.intermine.bio.util.PostProcessUtil;

import org.intermine.metadata.ClassDescriptor;
import org.intermine.metadata.ConstraintOp;
import org.intermine.metadata.MetaDataException;

import org.intermine.model.bio.Chromosome;
import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.DataSource;
import org.intermine.model.bio.Gene;
import org.intermine.model.bio.GeneFlankingRegion;
import org.intermine.model.bio.Location;
import org.intermine.model.bio.Supercontig;

import org.intermine.objectstore.ObjectStore;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.objectstore.ObjectStoreWriter;
import org.intermine.objectstore.query.ContainsConstraint;
import org.intermine.objectstore.query.Query;
import org.intermine.objectstore.query.QueryClass;
import org.intermine.objectstore.query.QueryCollectionReference;
import org.intermine.objectstore.query.Results;
import org.intermine.objectstore.query.ResultsRow;

import org.intermine.util.DynamicUtil;


import org.intermine.postprocess.PostProcessor;

/**
 * Create features to represent flanking regions of configurable distance either side of gene
 * features.  These will be used in overlap queries.
 *
 * @author rns
 * @author Sam Hokin
 */
public class CreateGeneFlankingFeaturesProcess extends PostProcessor {
    private static final Logger LOG = Logger.getLogger(CreateGeneFlankingFeaturesProcess.class);

    private ObjectStore os;
    private DataSet dataSet;
    private DataSource dataSource;

    // The sizes in kb of flanking regions to create.
    private static double[] distances = new double[] {0.5, 1, 2, 5, 10};

    // The values strings for up/down stream from a gene.
    private static String[] directions = new String[] {"upstream", "downstream"};

    Map<String,GeneFlankingRegion> flankingRegions = new ConcurrentHashMap<>(); // keyed by primaryIdentifier

    /**
     * Create a new instance
     *
     * @param osw object store writer
     */
    public CreateGeneFlankingFeaturesProcess(ObjectStoreWriter osw) throws ObjectStoreException {
        super(osw);
        this.os = osw.getObjectStore();
        dataSource = (DataSource) DynamicUtil.createObject(Collections.singleton(DataSource.class));
        dataSource.setName("InterMine post-processor");
        DataSource existingDataSource = os.getObjectByExample(dataSource, Collections.singleton("name"));
        if (existingDataSource==null) {
            // store new DataSource
            osw.beginTransaction();
            osw.store(dataSource);
            osw.commitTransaction();
        } else {
            // use existing DataSource
            dataSource = existingDataSource;
        }
        // DataSet - store later if we store any features
        dataSet = (DataSet) DynamicUtil.createObject(Collections.singleton(DataSet.class));
        dataSet.setName("InterMine gene-flanking regions");
        dataSet.setDescription("Gene-flanking regions created by the core InterMine post-processor");
        dataSet.setVersion("" + new Date()); // current time and date
        dataSet.setUrl("http://www.intermine.org");
        dataSet.setDataSource(dataSource);
    }

    /**
     * {@inheritDoc}
     */
    public void postProcess() throws ObjectStoreException {
        // query Genes that lack flanking regions
        Query q = new Query();
        QueryClass qcGene = new QueryClass(Gene.class);
        q.addFrom(qcGene);
        q.addToSelect(qcGene);
        QueryCollectionReference flankingRegionsRef = new QueryCollectionReference(qcGene, "flankingRegions");
        q.setConstraint(new ContainsConstraint(flankingRegionsRef, ConstraintOp.IS_NULL));
        Set<Gene> genes = new HashSet<>();
        Results results = os.execute(q, 1000, true, true, true);
        if (results.asList().size() > 0) {
            for (Object obj : results.asList()) {
                ResultsRow rr = (ResultsRow) obj;
                Gene gene = (Gene) rr.get(0);
                genes.add(gene);
            }
            LOG.info("Found " + genes.size() + " genes that need flanking regions.");
            /////////////////////////////////////////////
            // get flanking regions in a parallel stream
            genes.parallelStream().forEach(gene -> {
                    try {
                        createFlankingRegions(gene);
                    } catch (ObjectStoreException ex) {
                        System.err.println(ex);
                        System.exit(1);
                    }
                });
            /////////////////////////////////////////////
            LOG.info("Created " + flankingRegions.size() + " gene flanking regions.");
            LOG.info("Now storing...");
            // store
            osw.beginTransaction();
            osw.store(dataSet);
            for (GeneFlankingRegion flankingRegion : flankingRegions.values()) {
                Location location = null;
                if (flankingRegion.getChromosome() != null) {
                    location = flankingRegion.getChromosomeLocation();
                } else {
                    location = flankingRegion.getSupercontigLocation();
                }
                osw.store(flankingRegion);
                osw.store(location);
            }
            osw.commitTransaction();
            LOG.info("...done.");
        } else {
            LOG.info("All genes have flanking regions, none created.");
        }
    }

    /**
     * Create flanking regions for a gene, which are stored in the instance map.
     *
     * @param gene the Gene
     */
    void createFlankingRegions(Gene gene) throws ObjectStoreException {
        // get the location of the gene
        boolean onChromosome = (gene.getChromosome() != null);
        Chromosome chromosome = null;
        Supercontig supercontig = null;
        Location geneLoc = null;
        if (onChromosome) {
            chromosome = gene.getChromosome();
            geneLoc = gene.getChromosomeLocation();
        } else {
            supercontig = gene.getSupercontig();
            geneLoc = gene.getSupercontigLocation();
        }
        // run through the various flanking region distances
        for (double distance : distances) {
            // run through upstream and downstream
            for (String direction : directions) {
                String strand = geneLoc.getStrand();
                int geneStart = geneLoc.getStart().intValue();
                int geneEnd = geneLoc.getEnd().intValue();
                int contigLength = 0;
                if (onChromosome) {
                    contigLength = chromosome.getLength().intValue();
                } else {
                    contigLength = supercontig.getLength().intValue();
                }
                // gene touches a contig end so there isn't a flanking region
                if ((geneStart <= 1) || (geneEnd >= contigLength)) {
                    continue;
                }
                // create this GeneFlankingRegion
                GeneFlankingRegion region = (GeneFlankingRegion) DynamicUtil.createObject(Collections.singleton(GeneFlankingRegion.class));
                Location location = (Location) DynamicUtil.createObject(Collections.singleton(Location.class));
                region.setDistance(distance + "kb");
                region.setDirection(direction);
                try {
                    PostProcessUtil.checkFieldExists(os.getModel(), "GeneFlankingRegion", "includeGene", "Not setting");
                    region.setFieldValue("includeGene", false);
                } catch (MetaDataException e) {
                    // GeneFlankingRegion.includeGene not in model so do nothing
                }
                region.setGene(gene);
                if (onChromosome) {
                    region.setChromosome(chromosome);
                    region.setChromosomeLocation(location);
                } else {
                    region.setSupercontig(supercontig);
                    region.setSupercontigLocation(location);
                }
                region.setOrganism(gene.getOrganism());
                region.setStrain(gene.getStrain());
                String primaryIdentifier = gene.getPrimaryIdentifier() + "_" + distance + "_kb_" + direction;
                region.setPrimaryIdentifier(primaryIdentifier);
                // this should be some clever algorithm
                int start, end;
                if ("upstream".equals(direction) && "1".equals(strand)) {
                    start = geneStart - (int) Math.round(distance * 1000);
                    end = geneStart - 1;
                } else if ("upstream".equals(direction) && "-1".equals(strand)) {
                    start = geneEnd + 1;
                    end = geneEnd + (int) Math.round(distance * 1000);
                } else if ("downstream".equals(direction) && "1".equals(strand)) {
                    start = geneEnd + 1;
                    end = geneEnd + (int) Math.round(distance * 1000);
                } else {  // "downstream".equals(direction) && strand.equals("-1")
                    start = geneStart - (int) Math.round(distance * 1000);
                    end = geneStart - 1;
                }
                // if the region hangs off the start or end of a chromosome set it to finish
                // at the end of the chromosome
                location.setStart(Math.max(start, 1));
                int e = Math.min(end, contigLength);
                location.setEnd(e);
                location.setStrand(strand);
                if (onChromosome) {
                    location.setLocatedOn(chromosome);
                } else {
                    location.setLocatedOn(supercontig);
                }
                location.setFeature(region);
                region.setLength((location.getEnd().intValue() - location.getStart().intValue()) + 1);
                // store in our map
                flankingRegions.put(primaryIdentifier, region);
            }
        }
    }

}
