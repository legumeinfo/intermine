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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.intermine.bio.util.Constants;
import org.intermine.metadata.ConstraintOp;
import org.intermine.model.bio.OntologyAnnotation;
import org.intermine.model.bio.OntologyEvidence;
import org.intermine.model.bio.OntologyAnnotationEvidenceCode;
import org.intermine.model.bio.Gene;
import org.intermine.model.bio.OntologyTerm;
import org.intermine.model.bio.Protein;
import org.intermine.model.bio.Publication;
import org.intermine.bio.util.PostProcessUtil;
import org.intermine.objectstore.ObjectStore;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.objectstore.ObjectStoreWriter;
import org.intermine.objectstore.intermine.ObjectStoreInterMineImpl;
import org.intermine.objectstore.query.ConstraintSet;
import org.intermine.objectstore.query.ContainsConstraint;
import org.intermine.objectstore.query.Query;
import org.intermine.objectstore.query.QueryClass;
import org.intermine.objectstore.query.QueryCollectionReference;
import org.intermine.objectstore.query.QueryObjectReference;
import org.intermine.objectstore.query.Results;
import org.intermine.objectstore.query.ResultsRow;
import org.intermine.postprocess.PostProcessor;

/**
 * Take any OntologyAnnotation objects assigned to proteins and copy them to corresponding genes.
 * Merge evidence where duplication is found.
 * Update evidence codes with names and descriptions.
 * @author Richard Smith
 * @author julie sullivan
 */
public class GoPostprocess extends PostProcessor
{
    private static final Logger LOG = Logger.getLogger(GoPostprocess.class);
    protected ObjectStore os;

    /**
     * Create a new UpdateOrthologes object from an ObjectStoreWriter
     * @param osw writer on genomic ObjectStore
     */
    public GoPostprocess(ObjectStoreWriter osw) {
        super(osw);
        this.os = osw.getObjectStore();
    }

    /**
     * Copy all GO annotations from the Protein objects to the corresponding Gene(s)

     *
     * @throws ObjectStoreException if anything goes wrong
     */
    @Override
    public void postProcess() throws ObjectStoreException {

        long startTime = System.currentTimeMillis();

        osw.beginTransaction();

        Iterator<?> resIter = findProteinProperties(false);

        int count = 0;
        Gene lastGene = null;
        Map<OntologyTerm, OntologyAnnotation> annotations = new HashMap<OntologyTerm, OntologyAnnotation>();

        while (resIter.hasNext()) {
            ResultsRow<?> rr = (ResultsRow<?>) resIter.next();
            Gene thisGene = (Gene) rr.get(0);
            OntologyAnnotation thisAnnotation = (OntologyAnnotation) rr.get(1);

            // process last set of annotations if this is a new gene
            if (lastGene != null && !(lastGene.equals(thisGene))) {
                for (OntologyAnnotation item : annotations.values()) {
                    osw.store(item);
                }
                lastGene.setOntologyAnnotations(new HashSet(annotations.values()));
                LOG.debug("store gene " + lastGene.getSecondaryIdentifier() + " with "
                        + lastGene.getOntologyAnnotations().size() + " GO.");
                osw.store(lastGene);

                lastGene = thisGene;
                annotations = new HashMap<OntologyTerm, OntologyAnnotation>();
            }

            OntologyTerm term = thisAnnotation.getOntologyTerm();
            Set<OntologyEvidence> evidence = thisAnnotation.getEvidence();

            OntologyAnnotation tempAnnotation;
            try {
                tempAnnotation = PostProcessUtil.copyInterMineObject(thisAnnotation);
            } catch (IllegalAccessException e) {
                throw new RuntimeException(e);
            }

            if (hasDupes(annotations, term, evidence, tempAnnotation)) {
                // if a dupe, merge with already created object instead of creating new
                continue;
            }
            tempAnnotation.setSubject(thisGene);

            lastGene = thisGene;
            count++;
        }

        if (lastGene != null) {
            for (OntologyAnnotation item : annotations.values()) {
                osw.store(item);
            }
            lastGene.setOntologyAnnotations(new HashSet(annotations.values()));
            LOG.debug("store gene " + lastGene.getSecondaryIdentifier() + " with "
                    + lastGene.getOntologyAnnotations().size() + " GO.");
            osw.store(lastGene);
        }

        LOG.info("Created " + count + " new OntologyAnnotation objects for Genes"
                + " - took " + (System.currentTimeMillis() - startTime) + " ms.");
        osw.commitTransaction();
    }

    private boolean hasDupes(Map<OntologyTerm, OntologyAnnotation> annotations, OntologyTerm term,
            Set<OntologyEvidence> evidence, OntologyAnnotation newAnnotation) {
        boolean isDupe = false;
        OntologyAnnotation alreadySeenAnnotation = annotations.get(term);
        if (alreadySeenAnnotation != null) {
            isDupe = true;
            mergeEvidence(evidence, alreadySeenAnnotation);
        } else {
            annotations.put(term, newAnnotation);
        }
        return isDupe;
    }

    // we've seen this term, merge instead of storing new object
    private void mergeEvidence(Set<OntologyEvidence> evidence, OntologyAnnotation alreadySeenAnnotation) {
        for (OntologyEvidence g : evidence) {
            OntologyAnnotationEvidenceCode c = g.getCode();
            Set<Publication> pubs = g.getPublications();
            boolean foundMatch = false;
            for (OntologyEvidence alreadySeenEvidence : alreadySeenAnnotation.getEvidence()) {
                OntologyAnnotationEvidenceCode alreadySeenCode = alreadySeenEvidence.getCode();
                Set<Publication> alreadySeenPubs = alreadySeenEvidence.getPublications();
                // we've already seen this evidence code, just merge pubs
                if (c.equals(alreadySeenCode)) {
                    foundMatch = true;
                    alreadySeenPubs = mergePubs(alreadySeenPubs, pubs);
                }
            }
            if (!foundMatch) {
                // we don't have this evidence code
                alreadySeenAnnotation.addEvidence(g);
            }
        }
    }

    private Set<Publication> mergePubs(Set<Publication> alreadySeenPubs, Set<Publication> pubs) {
        Set<Publication> newPubs = new HashSet<Publication>();
        if (alreadySeenPubs != null) {
            newPubs.addAll(alreadySeenPubs);
        }
        if (pubs != null) {
            newPubs.addAll(pubs);
        }
        return newPubs;
    }

    /**
     * Query Gene->Protein->Annotation->OntologyTerm and return an iterator over the Gene,
     *  Protein and OntologyTerm.
     *
     * @param restrictToPrimaryGoTermsOnly Only get primary Annotation items linking the gene
     *  and the go term.
     */
    private Iterator<?> findProteinProperties(boolean restrictToPrimaryGoTermsOnly)
        throws ObjectStoreException {
        Query q = new Query();

        q.setDistinct(false);

        QueryClass qcGene = new QueryClass(Gene.class);
        q.addFrom(qcGene);
        q.addToSelect(qcGene);
        q.addToOrderBy(qcGene);

        QueryClass qcProtein = new QueryClass(Protein.class);
        q.addFrom(qcProtein);

        QueryClass qcAnnotation = new QueryClass(OntologyAnnotation.class);
        q.addFrom(qcAnnotation);
        q.addToSelect(qcAnnotation);

        ConstraintSet cs = new ConstraintSet(ConstraintOp.AND);

        QueryCollectionReference geneProtRef = new QueryCollectionReference(qcProtein, "genes");
        cs.addConstraint(new ContainsConstraint(geneProtRef, ConstraintOp.CONTAINS, qcGene));

        QueryObjectReference annSubjectRef =
            new QueryObjectReference(qcAnnotation, "subject");
        cs.addConstraint(new ContainsConstraint(annSubjectRef, ConstraintOp.CONTAINS, qcProtein));

        q.setConstraint(cs);

        ((ObjectStoreInterMineImpl) os).precompute(q, Constants.PRECOMPUTE_CATEGORY);
        Results res = os.execute(q, 5000, true, true, true);
        return res.iterator();
    }
}
