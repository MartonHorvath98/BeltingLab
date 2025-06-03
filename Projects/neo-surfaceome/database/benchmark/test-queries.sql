/* 1) Select a specific protein's info */
SELECT * FROM protein_info WHERE uniprot_id = 'P00533';

/* 2) List all domain names matching a keyword (e.g. kinase) */
SELECT * FROM domain_info WHERE domain_name LIKE '%kinase%';

/* 3) Count all transcripts */
SELECT COUNT(*) FROM gene_info;

/* 4) Get sequences longer than 1000 AA */
SELECT * FROM protein_seq WHERE LENGTH(sequence) > 1000;

/* 5) List SURFME genes */
SELECT * FROM surfme_filter;

/* 6) Top 10 longest transcripts */
SELECT transcript_id, _length FROM gene_info ORDER BY _length DESC LIMIT 10;

/* 7) Most common chromosome in gene_info */
SELECT chromosome, COUNT(*) FROM gene_info GROUP BY chromosome ORDER BY COUNT(*) DESC;

/* 8) Domain frequency */
SELECT domain_id, COUNT(*) AS domain_count FROM domain_annotation GROUP BY domain_id ORDER BY domain_count DESC LIMIT 10;

/* 9) Proteins with missing UniProt IDs */
SELECT protein_id FROM protein_info WHERE uniprot_id IS NULL;

/* 10) Proteins with the highest number of interaction partners */
SELECT protein_A, COUNT(*) FROM protein_interactions GROUP BY protein_A ORDER BY COUNT(*) DESC LIMIT 10;

/* 11) Get protein name and sequence */
SELECT pi.protein_name, ps.sequence
FROM protein_info pi
JOIN protein_seq ps ON pi.protein_id = ps.protein_id;

/* 12) List domains for a specific transcript */
SELECT da.transcript_id, da.domain_start, da.domain_end, di.domain_id, di.domain_name, di.annotation
FROM domain_annotation da
JOIN domain_info di ON da.domain_id = di.domain_id
WHERE da.transcript_id = 'ENST00000651593'
ORDER BY da.domain_start, da.domain_end DESC;

/* 13) Gene info of protein coding SURFME genes */
SELECT g.*
FROM gene_info g
JOIN surfme_filter s ON g.gene_id = s.gene_id
WHERE g.transcript_type  = 'protein_coding';

/* 14) Get protein interaction scores with names */
SELECT pi1.protein_name AS ProteinA, pi2.protein_name AS ProteinB, pi.combined_score
FROM protein_interactions pi
JOIN protein_info pi1 ON pi.protein_A = pi1.protein_id
JOIN protein_info pi2 ON pi.protein_B = pi2.protein_id
WHERE pi.combined_score > 900;

/* 15) Isoforms of genes with more than 2 domains */
SELECT gi.gene_symbol, da.transcript_id, COUNT(*)
FROM domain_annotation da
JOIN gene_info gi ON gi.transcript_id = da.transcript_id
GROUP BY da.transcript_id
HAVING COUNT(da.domain_id) > 2
ORDER BY gi.gene_symbol;

/* 16) Transcript length vs domain count */
SELECT gi.gene_symbol, gi.transcript_id, gi._length AS transcript_length, COUNT(da.domain_id) AS domain_count
FROM gene_info gi
JOIN domain_annotation da ON gi.transcript_id = da.transcript_id
GROUP BY gi.transcript_id
ORDER BY transcript_length DESC
LIMIT 10;

/* 17) IDs of surface proteins */
SELECT DISTINCT gi.protein_id
FROM gene_info gi
JOIN surfme_filter sf ON gi.gene_id = sf.gene_id
WHERE gi.protein_id != 'NA';

/* 18) Strong fusion partner proteins */
SELECT pi1.protein_name AS ProteinA, pi2.protein_name AS ProteinB, pi.combined_score
FROM protein_interactions pi
JOIN protein_info pi1 ON pi.protein_A = pi1.protein_id
JOIN protein_info pi2 ON pi.protein_B = pi2.protein_id
WHERE pi.combined_score >= 900 AND fusion_score = GREATEST(
		neighborhood_score,
		cooccurence_score,
		coexpression_score,
		experimental_score,
		database_score,
		textmining_score,
		fusion_score
	);

/* 19) Get significant domain annotation predictions, p-adj < 0.05 */
SELECT da.transcript_id, di.domain_name, da.score
FROM domain_annotation da
JOIN domain_info di ON da.domain_id = di.domain_id
WHERE da.accuracy >= 0.9 AND da.E_value< 0.05;

/* 20) Gene symbol and corresponding protein name */
SELECT gi.gene_symbol, pi.protein_name
FROM gene_info gi
JOIN protein_info pi ON gi.protein_id = pi.protein_id;

/* 21) Interaction partners of a protein, and their interacting domains */
SELECT DISTINCT gi1.gene_symbol AS GeneA, gi2.gene_symbol AS GeneB, di1.domain_id AS DomainA, di2.domain_id AS DomainB
FROM protein_interactions pi
JOIN gene_info gi1 ON pi.protein_A = gi1.protein_id
JOIN gene_info gi2 ON pi.protein_B = gi2.protein_id
JOIN domain_annotation di1 ON gi1.transcript_id = di1.transcript_id
JOIN domain_annotation di2 ON gi2.transcript_id = di2.transcript_id
WHERE pi.protein_A = 'ENSP00000269305' -- TP53
  AND di1.accuracy >= 0.9 AND di1.E_value < 0.05
  AND di2.accuracy >= 0.9 AND di2.E_value < 0.05
  AND pi.combined_score >= 900
  AND EXISTS (
      SELECT 1 FROM domain_interactions d
      WHERE d.domain_A = di1.domain_id AND d.domain_B = di2.domain_id
  )
LIMIT 10;

/* 22) All interactions involving surfme genes with domain-level details */
SELECT DISTINCT gi1.gene_symbol AS GeneA, gi2.gene_symbol AS GeneB, di1.domain_id AS DomainA, di2.domain_id AS DomainB
FROM protein_interactions pi
JOIN gene_info gi1 ON pi.protein_A = gi1.protein_id
JOIN gene_info gi2 ON pi.protein_B = gi2.protein_id
JOIN domain_annotation di1 ON gi1.transcript_id = di1.transcript_id
JOIN domain_annotation di2 ON gi2.transcript_id = di2.transcript_id
WHERE gi1.gene_id IN (
 		SELECT gene_id FROM surfme_filter WHERE gene_id IS NOT NULL
        )
	AND di1.accuracy >= 0.9 AND di1.E_value < 0.05
	AND di2.accuracy >= 0.9 AND di2.E_value < 0.05
	AND pi.combined_score >= 900
	AND EXISTS (
		SELECT 1 FROM domain_interactions d
		WHERE d.domain_A = di1.domain_id AND d.domain_B = di2.domain_id
	);

/* 23) Top 5 domains most often involved in interactions */
SELECT da.domain_id, d.annotation, COUNT(*) AS interaction_count
FROM domain_annotation da
JOIN domain_info d ON da.domain_id = d.domain_id
JOIN domain_interactions di ON da.domain_id = di.domain_A
WHERE da.accuracy >= 0.9 AND da.E_value < 0.05
GROUP BY da.domain_id
ORDER BY interaction_count DESC
LIMIT 10;

/* 24) How many interaction partners each domain of a protein has */
SELECT di1.domain_id AS DomainA, COUNT(*)
FROM protein_interactions pi
JOIN gene_info gi1 ON pi.protein_A = gi1.protein_id
JOIN gene_info gi2 ON pi.protein_B = gi2.protein_id
JOIN domain_annotation di1 ON gi1.transcript_id = di1.transcript_id
JOIN domain_annotation di2 ON gi2.transcript_id = di2.transcript_id
WHERE pi.protein_A = 'ENSP00000269305' -- TP53
	AND pi.combined_score >= 900
    AND di1.accuracy >= 0.9 AND di1.E_value < 0.05
    AND di2.accuracy >= 0.9 AND di2.E_value < 0.05
    AND EXISTS (
		SELECT 1 FROM domain_interactions d
		WHERE d.domain_A = di1.domain_id AND d.domain_B = di2.domain_id
	)
GROUP BY pi.protein_A, di1.domain_id
ORDER BY COUNT(*) DESC;

/* 25) All interacting domains for a given transcript */
SELECT da1.transcript_id AS source_transcript, di1.domain_name AS source_domain, di2.domain_name AS interacting_domain
FROM domain_annotation da1
JOIN domain_interactions d ON da1.domain_id = d.domain_A
JOIN domain_info di1 ON da1.domain_id = di1.domain_id
JOIN domain_info di2 ON d.domain_B = di2.domain_id
WHERE da1.transcript_id = 'ENST00000415913' -- IDH1
	AND  da1.accuracy >= 0.9 AND da1.E_value < 0.05; 

/* 26) Proteins from surfme genes and their interaction partners from non-surfme genes */
SELECT DISTINCT piA.protein_id AS surfme_protein, piB.protein_id AS other_protein
FROM gene_info giA
JOIN surfme_filter sf ON giA.gene_id = sf.gene_id
JOIN protein_interactions pi ON giA.protein_id = pi.protein_A
JOIN gene_info giB ON pi.protein_B = giB.protein_id
LEFT JOIN surfme_filter sf2 ON giB.gene_id = sf2.gene_id
JOIN protein_info piA ON giA.protein_id = piA.protein_id
JOIN protein_info piB ON giB.protein_id = piB.protein_id
WHERE sf2.gene_id IS NULL AND pi.combined_score >= 900;

/* 27) Full annotation pipeline for a transcript (gene + protein + domains) */
SELECT gi.gene_id, gi.gene_symbol, di.domain_id, di.domain_name, da.domain_start, da.domain_end
FROM gene_info gi
JOIN protein_info pi ON gi.protein_id = pi.protein_id
JOIN domain_annotation da ON gi.transcript_id = da.transcript_id
JOIN domain_info di ON da.domain_id = di.domain_id;


/* 28) Interaction partners of a protein on one of its domains */
SELECT DISTINCT gi1.gene_symbol AS GeneA, gi2.gene_symbol AS GeneB, di1.domain_id AS DomainA, di2.domain_id AS DomainB
FROM protein_interactions pi
JOIN gene_info gi1 ON pi.protein_A = gi1.protein_id
JOIN gene_info gi2 ON pi.protein_B = gi2.protein_id
JOIN domain_annotation di1 ON gi1.transcript_id = di1.transcript_id
JOIN domain_annotation di2 ON gi2.transcript_id = di2.transcript_id
WHERE pi.protein_A = 'ENSP00000269305' -- TP53
  AND di1.domain_id = 'PF00870' -- P53
  AND di1.accuracy >= 0.9 AND di1.E_value < 0.05
  AND di2.accuracy >= 0.9 AND di2.E_value < 0.05
  AND pi.combined_score >= 900
  AND EXISTS (
      SELECT 1 FROM domain_interactions d
      WHERE d.domain_A = di1.domain_id AND d.domain_B = di2.domain_id
  )
LIMIT 10;

/* 29) Domain–domain interactions potentially affected by a mutation */
SELECT gi1.gene_symbol, di1.domain_id AS affected_domain,
    di2.domain_id AS potentially_disrupted_ddi,
    gi2.gene_symbol AS potentially_disrupted_partner
FROM protein_interactions pi
JOIN gene_info gi1 ON pi.protein_A = gi1.protein_id
JOIN gene_info gi2 ON pi.protein_B = gi2.protein_id
JOIN domain_annotation di1 ON gi1.transcript_id = di1.transcript_id
JOIN domain_annotation di2 ON gi2.transcript_id = di2.transcript_id
WHERE gi1.transcript_id = 'ENST00000257290'  -- PDGFRA
	AND EXISTS (
		  SELECT 1 FROM domain_interactions d
		  WHERE d.domain_A = di1.domain_id AND d.domain_B = di2.domain_id
	)
	AND 150 BETWEEN di1.domain_start AND di1.domain_end; -- in PK_Tyr_Ser-Thr domain

/* 30) Overview of surfaceome proteins with domain and interaction stats */
SELECT 
    gi.gene_symbol,
    gi.transcript_id,
    gi.protein_id,
    ps._length AS protein_length,
    COUNT(DISTINCT da.domain_id) AS num_domains,
    COUNT(DISTINCT pi2.protein_B) AS num_interaction_partners
FROM surfme_filter sf
JOIN gene_info gi ON gi.gene_id = sf.gene_id
JOIN protein_seq ps ON gi.protein_id = ps.protein_id
LEFT JOIN domain_annotation da ON gi.transcript_id = da.transcript_id
LEFT JOIN protein_interactions pi2 ON gi.protein_id = pi2.protein_A
WHERE da.accuracy >= 0.9 AND da.E_Value < 0.05 AND pi2.combined_score >= 900
GROUP BY gi.gene_symbol, gi.transcript_id, gi.protein_id, ps._length
ORDER BY num_interaction_partners DESC, num_domains DESC
LIMIT 10;

