/* FREQUENT QUERY (1):
 -- Get all interacting domains for a given transcript --
 */
 -- 1) Pre-compute filtered, high-confidence domain annotation 
CREATE TEMPORARY TABLE tmp_high_conf_domain_annotation AS
SELECT * 
FROM domain_annotation
WHERE accuracy >= 0.9 AND E_value < 0.05;

DELIMITER $$
CREATE PROCEDURE MapDomains ( IN GeneId nvarchar(15))
BEGIN
SELECT DISTINCT gi.transcript_id, gi.gene_symbol, pi.protein_size, di.domain_name, da.domain_start, da.domain_end
FROM gene_info gi
JOIN protein_info pi ON gi.protein_id = pi.protein_id
JOIN tmp_high_conf_domain_annotation da ON gi.transcript_id = da.transcript_id
JOIN domain_info di ON da.domain_id = di.domain_id
WHERE gi.gene_id = GeneId;
END$$
DELIMITER ;

/* Gene Table:                                                                                                                                                                          |
| ----------- | --------------- | ------------------ |
| Gene Symbol | Ensembl Gene ID | UniProt Protein ID |                                                                                                                                                                          |
| ----------- | --------------- | ------------------ |
| EGFR        | ENSG00000146648 | P00533             |                                                                                                                                                                          |
| PDPN        | ENSG00000162493 | Q86YL7             |                                                                                                                                                                          |
| CD44        | ENSG00000026508 | P16070             |                                                                                                                                                                          |
| ITGA6       | ENSG00000115221 | P23229             |                                                                                                                                                                          |
| L1CAM       | ENSG00000198910 | P32004             |                                                                                                                                                                          |
| PROM1       | ENSG00000007062 | O43490             |                                                                                                                                                                          |
| IL13RA2     | ENSG00000123496 | Q14627             |                                                                                                                                                                          |
| MET         | ENSG00000105976 | P08581             |                                                                                                                                                                          |
| CD276       | ENSG00000103855 | Q5ZPR3             |                                                                                                                                                                          |
| ERBB2       | ENSG00000141736 | P04626             |                                                                                                                                                                          |
| CD70        | ENSG00000125726 | P32970             |                                                                                                                                                                          |
| CEACAM1     | ENSG00000086548 | P13688             |                                                                                                                                                                          |
| PVR         | ENSG00000073008 | P15151             |                                                                                                                                                                          |
| FAP         | ENSG00000078098 | Q12884             |                                                                                                                                                                          |
| VCAM1       | ENSG00000162692 | P19320             |                                                                                                                                                                          |
| CLDN1       | ENSG00000163347 | O95832             |                                                                                                                                                                          |
| MUC1        | ENSG00000185499 | P15941             |                                                                                                                                                                          |
| TSPAN1      | ENSG00000165629 | O60635             |                                                                                                                                                                          |
| NTRK2       | ENSG00000148053 | Q16620             |                                                                                                                                                                          |
| ----------- | --------------- | ------------------ | */
CALL MapDomains("ENSG00000146648");
/* EGFR (1210aa)
Furin-like	10	149
Recep_L_domain	1	111
PK_Tyr_Ser-Thr	3	256
GF_recep_IV	1	131
*/
CALL MapDomains("ENSG00000162493");
/* PDPN (238aa)
Podoplanin	1	138
*/
CALL MapDomains("ENSG00000026508");
/* CD44 (742aa)
Xlink	2	93
*/
CALL MapDomains("ENSG00000115221");
/* ITGA6 (788aa)
Integrin_beta	2	248
Integrin_B_tail	1	85
EGF_2	1	32
Integrin_b_cyt	1	45
PSI_integrin	6	46
I-EGF_1	1	29 */
CALL MapDomains("ENSG00000198910");
/* L1CAM (1257aa)
ig	6	81
I-set	10	90
Bravo_FIGEY	1	88
*/
CALL MapDomains("ENSG00000007062");
/* PROM1 (865aa)
RasGAP_C	36	85
Prominin	2	799
*/
CALL MapDomains("ENSG00000123496");
/* IL13RA2 (380aa)
IL6Ra-bind	1	96
*/
CALL MapDomains("ENSG00000105976");
/* MET (1408aa)
TIG	1	85
PK_Tyr_Ser-Thr	2	258
*/
CALL MapDomains("ENSG00000103855");
/* CD276 - no hits */
CALL MapDomains("ENSG00000141736");
/* ERBB2 (1255aa)
PK_Tyr_Ser-Thr	3	258
Recep_L_domain	1	111
GF_recep_IV	2	131
Furin-like	10	149
*/
CALL MapDomains("ENSG00000125726");
/* CD70 - no hits */
CALL MapDomains("ENSG00000086548");
/* CEACAM6 (344aa)
Ig_2	2	76
V-set	3	108
ig	8	80
*/
CALL MapDomains("ENSG00000073008");
/* PVR (417aa)
C2-set_2	2	88
ig	12	84
*/
CALL MapDomains("ENSG00000078098");
/* FAP (760aa)
DPPIV_N	1	354
Peptidase_S9	3	210
DPPIV_rep	1	21
*/
CALL MapDomains("ENSG00000162692");
/* VCAM1 (739aa)
I-set	4	90
C2-set	1	80
ig	2	85
*/
CALL MapDomains("ENSG00000163347");
/* CLDN1 (211aa)
PMP22_Claudin	2	166
*/
CALL MapDomains("ENSG00000165629");
/* TSPAN1 (241aa)
Tetraspanin	1	229
*/
CALL MapDomains("ENSG00000148053");
/* NTRK2 (838aa)
PK_Tyr_Ser-Thr	2	258
TPKR_C2	2	45
ig	2	86
LRR_8	2	60
LRRNT	1	28
*/