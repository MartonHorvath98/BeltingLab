/* FREQUENT QUERY (2):
 -- Get interacting proteins affected by a mutation --
 */
 -- 1) Pre-compute filtered, high-confidence protein interactions 
CREATE TEMPORARY TABLE tmp_high_conf_protein_interactions AS
SELECT * 
FROM protein_interactions
WHERE combined_score >= 900;
-- 2) Precompute protein mapping to every protein
CREATE TABLE protein_domain_map AS
SELECT gi.gene_id, pi.protein_id, di.domain_id, di.domain_name, da.domain_start, da.domain_end, da.accuracy, da.E_value
FROM gene_info gi
JOIN protein_info pi ON gi.protein_id = pi.protein_id
JOIN tmp_high_conf_domain_annotation da ON gi.transcript_id = da.transcript_id
JOIN domain_info di ON da.domain_id = di.domain_id;
-- 3) Materialize protein-domain-domain-protein interaction map
CREATE TABLE domain_domain_protein_map AS
SELECT DISTINCT pdm1.gene_id AS GeneA, pdm1.domain_id AS DomainA, pdm2.domain_id AS DomainB, pdm2.gene_id AS GeneB 
FROM tmp_high_conf_protein_interactions pi
JOIN protein_domain_map pdm1 ON pi.protein_A = pdm1.protein_id
JOIN protein_domain_map pdm2 ON pi.protein_B = pdm2.protein_id
JOIN domain_interactions ddi 
  ON ddi.domain_A = pdm1.domain_id AND ddi.domain_B = pdm2.domain_id;
  
DELIMITER $$
CREATE PROCEDURE PerturbedDomains ( IN GeneId nvarchar(15), SnvLoc int)
BEGIN
SELECT DISTINCT gi1.gene_symbol AS affected_gene,
	ddi.DomainA AS affected_domain,
    ddi.DomainB AS potentially_disrupted_ddi,
    gi2.gene_symbol AS potentially_disrupted_partner
FROM protein_domain_map pdm
JOIN domain_domain_protein_map ddi ON pdm.gene_id = ddi.GeneA AND pdm.domain_id = ddi.DomainA
JOIN gene_info gi1 ON ddi.GeneA = gi1.gene_id
JOIN gene_info gi2 ON ddi.GeneB = gi2.gene_id
WHERE pdm.gene_id = GeneId
	AND SnvLoc BETWEEN pdm.domain_start AND pdm.domain_end;
END$$
DELIMITER ;
-- DROP procedure PerturbedDomains;

CALL PerturbedDomains("ENSG00000146648",150);
-- 1) EGFR (150) - 95 row(s) returned 0.313 sec [e.g. PTK2, PIK3R2, JAK2, etc.]
CALL PerturbedDomains("ENSG00000162493",150);
-- 2) PDPN (70) - 0 row(s) returned 0.219 sec
CALL PerturbedDomains("ENSG00000026508",70);
-- 3) CD44 (70) - 2 row(s) returned 0.641 sec [e.g. TNFAIP6, VCAN]
CALL PerturbedDomains("ENSG00000115221",200);
-- 4) ITGA6 (200) - 69 row(s) returned 0.297 sec [e.g. TNR, FN1, ITGB2, etc.]
CALL PerturbedDomains("ENSG00000198910",89);
-- 5) L1CAM (89) - 19 row(s) returned 0.375 sec [e.g. ANK1, NCAM2, CNTN1, etc.]
CALL PerturbedDomains("ENSG00000007062",200);
-- 6) PROM1 (200) - 0 row(s) returned 0.422
CALL PerturbedDomains("ENSG00000123496",70);
-- 7) IL13RA2 (70) - 5 row(s) returned 0.187 sec [e.g. IL4, IL13, IL13RA1]
CALL PerturbedDomains("ENSG00000105976",70);
-- 8) MET (70) - 49 row(s) returned 0.250 sec [e.g. SHC1, PTK2, JAK2]
CALL PerturbedDomains("ENSG00000141736",200);
-- 9) ERBB2 (200) - 42 row(s) returned 0.484 sec [e.g. PTK2, SHC1, JAK2]
CALL PerturbedDomains("ENSG00000086548",90);
-- 10) CEACAM6 (90) - 3 row(s) returned 0.172 sec [e.g. CEACAM8]
CALL PerturbedDomains("ENSG00000073008",70);
-- 11) PVR (70) - 4 row(s) returned 0.219 sec [e.g. CD96, NECTIN3 ]
CALL PerturbedDomains("ENSG00000078098",300);
-- 12) FAP (300) - 0 row(s) returned 0.328 sec
CALL PerturbedDomains("ENSG00000162692",40);
-- 13) VCAM1 (40) - 2 row(s) returned 0.218 sec [ FN1 ]
CALL PerturbedDomains("ENSG00000163347",150);
-- 14) CLDN1 (150) - 1 row(s) returned 0.157 sec [ CLDN4 ]
CALL PerturbedDomains("ENSG00000165629",150);
-- 15) TSPAN1 (150) - 0 row(s) returned 0.219 sec
CALL PerturbedDomains("ENSG00000148053",150);
-- 16) NTRK2 (150) - 29 row(s) returned 0.281 sec [e.g. SHC1, SHC3, PLCG2, etc.]