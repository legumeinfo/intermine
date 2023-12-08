package org.intermine.bio.postprocess;

/**
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.log4j.Logger;
import org.intermine.bio.util.ClobAccessReverseComplement;
import org.intermine.bio.util.Constants;
import org.intermine.metadata.ConstraintOp;
import org.intermine.bio.util.PostProcessUtil;
import org.intermine.metadata.MetaDataException;
import org.intermine.metadata.Model;
import org.intermine.model.bio.CDS;
import org.intermine.model.bio.Chromosome;
import org.intermine.model.bio.Exon;
import org.intermine.model.bio.Gene;
import org.intermine.model.bio.Location;
import org.intermine.model.bio.Sequence;
import org.intermine.model.bio.SequenceFeature;
import org.intermine.model.bio.Supercontig;
import org.intermine.model.bio.Transcript;
import org.intermine.objectstore.ObjectStore;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.objectstore.ObjectStoreWriter;
import org.intermine.objectstore.intermine.ObjectStoreInterMineImpl;
import org.intermine.objectstore.proxy.ProxyReference;
import org.intermine.objectstore.query.ClobAccess;
import org.intermine.objectstore.query.ClassConstraint;
import org.intermine.objectstore.query.ConstraintSet;
import org.intermine.objectstore.query.ContainsConstraint;
import org.intermine.objectstore.query.PendingClob;
import org.intermine.objectstore.query.Query;
import org.intermine.objectstore.query.QueryClass;
import org.intermine.objectstore.query.QueryCollectionReference;
import org.intermine.objectstore.query.QueryField;
import org.intermine.objectstore.query.QueryNode;
import org.intermine.objectstore.query.QueryObjectReference;
import org.intermine.objectstore.query.QueryValue;
import org.intermine.objectstore.query.Results;
import org.intermine.objectstore.query.ResultsRow;
import org.intermine.objectstore.query.SimpleConstraint;
import org.intermine.objectstore.query.SingletonResults;
import org.intermine.postprocess.PostProcessor;
import org.intermine.util.DynamicUtil;

/**
 * Transfer sequences from the contigs (Chromosome and Supercontig) to the SequenceFeature objects
 * which lack sequence.
 *
 * @author Kim Rutherford
 * @author Sam Hokin
 */
public class TransferSequencesProcess extends PostProcessor {

    private static final Logger LOG = Logger.getLogger(TransferSequencesProcess.class);

    private Model model; // stored for various methods

    /**
     * @param osw object store writer
     */
    public TransferSequencesProcess(ObjectStoreWriter osw) {
        super(osw);
    }

    /**
     * {@inheritDoc}
     */
    public void postProcess() throws IllegalAccessException, MetaDataException, ObjectStoreException {
        model = Model.getInstanceByName("genomic");
        transferToChromosomeFeatures();
        transferToSupercontigFeatures();
        transferToCDSes();
        transferToTranscripts();
    }

    /**
     * Use the Location relations to copy the sequence from the Chromosomes to the
     * SequenceFeatures located on them and which don't already have a sequence.
     * This method queries chromosomes, then populate features on them one by one.
     */
    protected void transferToChromosomeFeatures() throws IllegalAccessException, ObjectStoreException {
        long startTime = System.currentTimeMillis();
        // query Chromosomes
        Query q = new Query();
        QueryClass qc = new QueryClass(Chromosome.class);
        q.addFrom(qc);
        q.addToSelect(qc);
        QueryObjectReference seqRef = new QueryObjectReference(qc, "sequence");
        ContainsConstraint cc = new ContainsConstraint(seqRef, ConstraintOp.IS_NOT_NULL);
        q.setConstraint(cc);
        Set<Chromosome> chromosomes = new HashSet<>();
        SingletonResults results = osw.getObjectStore().executeSingleton(q);
        for (Object obj : results.asList()) {
            Chromosome chr = (Chromosome) obj;
            chromosomes.add(chr);
        }
        LOG.info("Found " + chromosomes.size() + " chromosomes with sequence, took " + (System.currentTimeMillis() - startTime) + " ms.");
        // do the transfer work for Chromosomes
        for (Chromosome chromosome : chromosomes) {
            int numFeatures = transferForChromosome(chromosome);
        }
    }

    /**
     * Use the Location relations to copy the sequence from Supercontigs to the
     * SequenceFeatures located on them and which don't already have a sequence.
     * This method queries features on supercontigs and populates them en masse since a given
     * Supercontig typically has zero to a few features.
     */
    protected void transferToSupercontigFeatures() throws IllegalAccessException, ObjectStoreException {
        long startTime = System.currentTimeMillis();
        LOG.info("Starting sequence transfer to features on supercontigs...");
        // query features on supercontigs
        Query q = new Query();
        ConstraintSet cs = new ConstraintSet(ConstraintOp.AND);
        QueryClass qcFeature = new QueryClass(SequenceFeature.class);
        q.addFrom(qcFeature);
        q.addToSelect(qcFeature);
        // only features with null sequence
        QueryObjectReference seqRef = new QueryObjectReference(qcFeature, "sequence");
        ContainsConstraint seqNullConstraint = new ContainsConstraint(seqRef, ConstraintOp.IS_NULL);
        cs.addConstraint(seqNullConstraint);
        // only features on Supercontigs
        QueryObjectReference scRef = new QueryObjectReference(qcFeature, "supercontig");
        ContainsConstraint scNotNullConstraint = new ContainsConstraint(scRef, ConstraintOp.IS_NOT_NULL);
        cs.addConstraint(scNotNullConstraint);
        q.setConstraint(cs);
        // run the query, storing the feature sequences as we go
        int count = 0;
        SingletonResults results = osw.getObjectStore().executeSingleton(q);
        osw.beginTransaction();
        for (Object obj : results.asList()) {
            SequenceFeature feature = (SequenceFeature) obj;
            if (PostProcessUtil.isInstance(model, feature, "CDS")) continue; // TO DO?
            if (PostProcessUtil.isInstance(model, feature, "Transcript")) continue; // TO DO?
            Supercontig supercontig = feature.getSupercontig();
            Location location = feature.getSupercontigLocation();
            // extract the feature's sequence from the supercontig
            ClobAccess featureSeq = getSubSequence(supercontig.getSequence(), location);
            if (featureSeq == null) {
                LOG.warn("Could not get feature sequence for location: " + location);
                continue;
            }
            // store the feature sequence and the feature clone
            Sequence sequence = (Sequence) DynamicUtil.createObject(Collections.singleton(Sequence.class));
            sequence.setResidues(featureSeq);
            sequence.setLength(featureSeq.length());
            osw.store(sequence);
            SequenceFeature clone = PostProcessUtil.cloneInterMineObject(feature);
            clone.setSequence(sequence);
            clone.setLength(featureSeq.length());
            osw.store(clone);
            count++;
            if (count % 10000 == 0) LOG.info(count + " sequences stored for supercontig features...");
        }
        osw.commitTransaction();
        LOG.info("Stored " + count + " sequences for features on supercontigs; took " + (System.currentTimeMillis() - startTime) + " ms.");
    }
    
    /**
     * Transfer sequences to features on the given Chromosome.
     *
     * @param chromosome
     * @return the number of features on this chromosome that lacked sequences
     */
    protected int transferForChromosome(Chromosome chromosome) throws IllegalAccessException, ObjectStoreException {
        long startTime = System.currentTimeMillis();

        // some constants
        int chromosomeId = chromosome.getId();
        String chromosomeIdentifier = chromosome.getPrimaryIdentifier();
        Sequence chromosomeSequence = chromosome.getSequence();

        // our query
        Query q = new Query();
        q.setDistinct(false);
        // Chromosome NOT SELECTED
        QueryClass qcChromosome = new QueryClass(Chromosome.class);
        q.addFrom(qcChromosome);
        // 0 SequenceFeature
        QueryClass qcFeature = new QueryClass(SequenceFeature.class);
        q.addFrom(qcFeature);
        q.addToSelect(qcFeature);
        // 1 Location
        QueryClass qcLocation = new QueryClass(Location.class);
        q.addFrom(qcLocation);
        q.addToSelect(qcLocation);

        // we have a boatload of AND constraints
        ConstraintSet cs = new ConstraintSet(ConstraintOp.AND);

        // Chromosome must have given chromosome id
        QueryField qfId = new QueryField(qcChromosome, "id");
        cs.addConstraint(new SimpleConstraint(qfId, ConstraintOp.EQUALS, new QueryValue(chromosomeId)));
        // Location must be located on our Chromosome
        QueryObjectReference ref1 = new QueryObjectReference(qcLocation, "locatedOn");
        ContainsConstraint cc1 = new ContainsConstraint(ref1, ConstraintOp.CONTAINS, qcChromosome);
        cs.addConstraint(cc1);
        // Location must be for the SequenceFeature
        QueryObjectReference ref2 = new QueryObjectReference(qcLocation, "feature");
        ContainsConstraint cc2 = new ContainsConstraint(ref2, ConstraintOp.CONTAINS, qcFeature);
        cs.addConstraint(cc2);
        // the feature sequence must be null
        QueryObjectReference lsfSeqRef = new QueryObjectReference(qcFeature, "sequence");
        ContainsConstraint lsfSeqRefNull = new ContainsConstraint(lsfSeqRef, ConstraintOp.IS_NULL);
        cs.addConstraint(lsfSeqRefNull);

        // set the constraint
        q.setConstraint(cs);

        // create precompute indexes on our QueryClass objects to speed things up
        Set<QueryNode> indexesToCreate = new HashSet<>();
        indexesToCreate.add(qcLocation);
        indexesToCreate.add(qcFeature);
        ((ObjectStoreInterMineImpl) osw.getObjectStore()).precompute(q, indexesToCreate, Constants.PRECOMPUTE_CATEGORY);

        // run the query
        Results results = osw.getObjectStore().execute(q, 1000, true, true, true);
        int numResults = results.size();
        if (numResults > 0) {
            int count = 0;
            osw.beginTransaction();
            for (Object obj : results.asList()) {
                ResultsRow rr = (ResultsRow) obj;
                SequenceFeature feature = (SequenceFeature) rr.get(0);
                Location location = (Location) rr.get(1);
                // these are done in transfertoCDSes and transferToTranscripts
                if (PostProcessUtil.isInstance(model, feature, "CDS")) continue; // done in transferToCDSes;
                if (PostProcessUtil.isInstance(model, feature, "Transcript")) continue; // done in transferToTranscripts;
                // bail on certain types of feature
                if (PostProcessUtil.isInstance(model, feature, "GeneticMarker")) continue;
                // bail if Gene too long, boss!
                if (feature instanceof Gene) {
                    Gene gene = (Gene) feature;
                    if (gene.getLength() != null && gene.getLength().intValue() > 2000000) {
                        LOG.warn("Gene too long to transfer sequence, ignoring: " + gene);
                        continue;
                    }
                }
                // get the feature's sequence
                ClobAccess featureSeq = getSubSequence(chromosomeSequence, location);
                if (featureSeq == null) {
                    // probably the location is out of range
                    LOG.warn("Could not get feature sequence for location: " + location);
                    continue;
                }
                // store the feature sequence and the feature clone
                Sequence sequence = (Sequence) DynamicUtil.createObject(Collections.singleton(Sequence.class));
                sequence.setResidues(featureSeq);
                sequence.setLength(featureSeq.length());
                osw.store(sequence);
                SequenceFeature clone = PostProcessUtil.cloneInterMineObject(feature);
                clone.setSequence(sequence);
                clone.setLength(featureSeq.length());
                osw.store(clone);
                count++;
            }
            osw.commitTransaction();
            LOG.info("Stored " + count + " feature sequences for " + chromosomeIdentifier + "; took " + (System.currentTimeMillis() - startTime) + " ms.");
        }
        // return full numResults since we skipped CDSes and transcripts
        return numResults;
    }

    /**
     * Get the subsequence for a Location on a given Sequence.
     *
     * @param contigSequence the Sequence
     * @param location the Location
     * @return the ClobAccess holding the subsequence spanned by the location
     */
    private static ClobAccess getSubSequence(Sequence contigSequence, Location location) {
        int charsToCopy = location.getEnd().intValue() - location.getStart().intValue() + 1;
        ClobAccess contigSequenceString = contigSequence.getResidues();
        if (charsToCopy > contigSequenceString.length()) {
            LOG.warn("SequenceFeature too long, ignoring; Location:" + location.getId() + " feature:" + location.getFeature());
            return null;
        }
        int startPos = location.getStart().intValue() - 1;
        int endPos = startPos + charsToCopy;
        if (startPos < 0 || endPos < 0) {
            LOG.warn("SequenceFeature has negative coordinate, ignoring Location:" + location.getId() + " feature:" + location.getFeature());
            return null;
        }
        if (endPos > contigSequenceString.length()) {
            LOG.warn("End coordinate is greater than chromsome length, ignoring Location:" + location.getId() + " feature:" + location.getFeature());
            return null;
        }
        ClobAccess subSeqString;
        if (startPos < endPos) {
            subSeqString = contigSequenceString.subSequence(startPos, endPos);
        } else {
            subSeqString = contigSequenceString.subSequence(endPos, startPos);
        }
        if ("-1".equals(location.getStrand())) {
            subSeqString = new ClobAccessReverseComplement(subSeqString);
        }
        return subSeqString;
    }

    /**
     * For every Transcript, join and transfer the sequences from the child Exons to a new Sequence
     * object for the transcript.
     */
    protected void transferToTranscripts() throws MetaDataException, ObjectStoreException {
        String message = "Transferring exon sequences to transcripts;";
        PostProcessUtil.checkFieldExists(model, "Transcript", "exons", message);
        PostProcessUtil.checkFieldExists(model, "Exon", null, message);

        long startTime = System.currentTimeMillis();

        // our query
        Query q = new Query();
        q.setDistinct(false);
        
        // 0 Transcript
        QueryClass qcTranscript = new QueryClass(Transcript.class);
        q.addFrom(qcTranscript);
        q.addToSelect(qcTranscript);

        // 1 Exon
        QueryClass qcExon = new QueryClass(Exon.class);
        q.addFrom(qcExon);
        q.addToSelect(qcExon);

        ConstraintSet cs = new ConstraintSet(ConstraintOp.AND);

        // Transcript has exons collection
        QueryCollectionReference exonsRef = new QueryCollectionReference(qcTranscript, "exons");
        ContainsConstraint transcriptHasExons = new ContainsConstraint(exonsRef, ConstraintOp.CONTAINS, qcExon);
        cs.addConstraint(transcriptHasExons);
        
        // Exon has sequence reference
        QueryObjectReference exonSeqRef = new QueryObjectReference(qcExon, "sequence");
        ContainsConstraint exonHasSequence = new ContainsConstraint(exonSeqRef, ConstraintOp.IS_NOT_NULL);
        cs.addConstraint(exonHasSequence);
        
        // Transcript does NOT have sequence
        QueryObjectReference transcriptSeqRef = new QueryObjectReference(qcTranscript, "sequence");
        ContainsConstraint transcriptLacksSequence = new ContainsConstraint(transcriptSeqRef, ConstraintOp.IS_NULL);
        cs.addConstraint(transcriptLacksSequence);

        q.setConstraint(cs);

        // create precompute indexes on our QueryClass objects to speed things up a shade
        Set<QueryNode> indexesToCreate = new HashSet<>();
        indexesToCreate.add(qcTranscript);
        indexesToCreate.add(qcExon);
        ((ObjectStoreInterMineImpl) osw.getObjectStore()).precompute(q, indexesToCreate, Constants.PRECOMPUTE_CATEGORY);

        Results results = osw.getObjectStore().execute(q, 1000, true, true, true);
        if (results.size() > 0) {
            int count = 0;
            int numChromosomeExons = 0;
            int numSupercontigExons = 0;

            // load transcripts into a map keyed by id
            // load subsequences into a map of TreeMaps keyed by id, with TreeMaps keyed by start
            Map<Integer,Transcript> transcripts = new HashMap<>();
            Map<Integer,TreeMap<Integer,String>> transcriptSequences = new HashMap<>();
            for (Object obj : results.asList()) {
                ResultsRow rr = (ResultsRow) obj;
                Transcript transcript =  (Transcript) rr.get(0);
                Exon exon = (Exon) rr.get(1);
                Location location = null;
                if (exon.getSupercontig() != null) {
                    location = exon.getSupercontigLocation();
                    numSupercontigExons++;
                } else {
                    location = exon.getChromosomeLocation();
                    numChromosomeExons++;
                }
                if (transcripts.containsKey(transcript.getId())) {
                    // add exon sequence to exon sequence TreeMap
                    TreeMap<Integer,String> exonSequences = transcriptSequences.get(transcript.getId());
                    exonSequences.put(location.getStart(), exon.getSequence().getResidues().toString());
                } else {
                    // add new transcript and initialize new exon sequence TreeMap
                    transcripts.put(transcript.getId(), transcript);
                    TreeMap<Integer,String> exonSequences = new TreeMap<>();
                    exonSequences.put(location.getStart(), exon.getSequence().getResidues().toString());
                    transcriptSequences.put(transcript.getId(), exonSequences);
                }
                count++;
            }
            LOG.info("Retrieved " + count + " Transcript exon sequences, " + numChromosomeExons + " on chromosomes, " + numSupercontigExons + " on supercontigs; " +
                     "took " + (System.currentTimeMillis() - startTime) + " ms.");

            // now splice together the sorted exon sequences for each transcript
            startTime = System.currentTimeMillis();
            count = 0;
            osw.beginTransaction();
            for (Integer id : transcripts.keySet()) {
                Transcript transcript = transcripts.get(id);
                TreeMap<Integer,String> exonSequences = transcriptSequences.get(id);
                StringBuffer transcriptSequence = new StringBuffer();
                for (String sequence : exonSequences.values()) {
                    transcriptSequence.append(sequence);
                }
                storeNewSequence(transcript, new PendingClob(transcriptSequence.toString()));
                count++;
            }
            osw.commitTransaction();
            LOG.info("Stored " + count + " Transcript sequences; took " + (System.currentTimeMillis() - startTime) + " ms.");
        } else {
            LOG.info("All transcript objects already contain sequence.");
        }
    }

    /**
     * For each CDS, join and transfer the sequences from each of the CDS.locations to a new Sequence
     * object stored for the CDS.
     */
    private void transferToCDSes() throws MetaDataException, ObjectStoreException {
        long startTime = System.currentTimeMillis();

        Query q = new Query();
        q.setDistinct(false);

        // CDS
        QueryClass qcCDS = new QueryClass(CDS.class);
        q.addFrom(qcCDS);
        q.addToSelect(qcCDS);
        q.addToOrderBy(qcCDS);

        // CDS does NOT have sequence
        QueryObjectReference sequenceRef = new QueryObjectReference(qcCDS, "sequence");
        ContainsConstraint sequenceIsNull = new ContainsConstraint(sequenceRef, ConstraintOp.IS_NULL);
        q.setConstraint(sequenceIsNull);

        // create precompute indexes on our QueryClass objects to speed things up
        Set<QueryNode> indexesToCreate = new HashSet<>();
        indexesToCreate.add(qcCDS);
        ((ObjectStoreInterMineImpl) osw.getObjectStore()).precompute(q, indexesToCreate, Constants.PRECOMPUTE_CATEGORY);
        
        // execute the query, stitching locations for each CDS and storing as we go
        SingletonResults results = osw.getObjectStore().executeSingleton(q);
        if (results.size() > 0) {
            osw.beginTransaction();
            int count = 0;
            int numChromosomeLocations = 0;
            int numSupercontigLocations = 0;
            Map<Integer,TreeMap<Integer,String>> cdsSequences = new HashMap<>();
            for (Object obj : results.asList()) {
                CDS cds = (CDS) obj;
                // put locations into start-sorted TreeMap
                TreeMap<Integer,String> sequences = new TreeMap<>();
                for (Location location : cds.getLocations()) {
                    if (cds.getSupercontig() != null) {
                        Supercontig contig = cds.getSupercontig();
                        String sequence = getSubSequence(contig.getSequence(), location).toString();
                        sequences.put(location.getStart(), sequence);
                        numSupercontigLocations++;
                    } else {
                        Chromosome contig = cds.getChromosome();
                        String sequence = getSubSequence(contig.getSequence(), location).toString();
                        sequences.put(location.getStart(), sequence);
                        numChromosomeLocations++;
                    }
                }
                // stich together the start-sorted sequences and store
                StringBuffer cdsSequence = new StringBuffer();
                for (String sequence : sequences.values()) {
                    cdsSequence.append(sequence);
                }
                storeNewSequence(cds, new PendingClob(cdsSequence.toString()));
                count++;
            }
            osw.commitTransaction();
            LOG.info("Stored " + count + " CDS sequences from " +
                     numChromosomeLocations + " chromosome and " + numSupercontigLocations + " supercontig locations; " +
                     "took " + (System.currentTimeMillis() - startTime) + " ms.");
        } else {
            LOG.info("All CDS objects already contain sequence.");
        }            
    }

    /**
     * Store a new Sequence along with its updated SequenceFeature.
     * Utility method used by transferToCDSes and transferToTranscripts.
     */
    private void storeNewSequence(SequenceFeature feature, ClobAccess sequenceString) throws ObjectStoreException {
        Sequence sequence = (Sequence) DynamicUtil.createObject(Collections.singleton(Sequence.class));
        sequence.setResidues(sequenceString);
        sequence.setLength(sequenceString.length());
        osw.store(sequence);
        feature.proxySequence(new ProxyReference(osw.getObjectStore(), sequence.getId(), Sequence.class));
        feature.setLength(sequenceString.length());
        osw.store(feature);
    }
    
}
