-- SET GLOBAL local_infile=1;
-- DROP TABLE IF EXISTS surfme_filter;

CREATE TABLE surfme_filter (
	gene_id VARCHAR(15) PRIMARY KEY
);

-- LOAD DATA LOCAL INFILE "C:\\Users\\gotar\\GitHub\\BeltingLab\\Projects\\neo-surfaceome\\database\\surfme_genes.txt"
LOAD DATA INFILE '/docker-entrypoint-initdb.d/data/surfme_genes.txt'
INTO TABLE surfme_filter
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(gene_id);
 
-- DROP TABLE IF EXISTS protein_info;
CREATE TABLE protein_info (
	protein_id VARCHAR(15) PRIMARY KEY, 
    protein_name VARCHAR(255),
    protein_size INT,
    annotation VARCHAR(255)
);

-- LOAD DATA LOCAL INFILE "C:\\Users\\gotar\\GitHub\\BeltingLab\\Projects\\neo-surfaceome\\database\\protein_info_clean.txt"
LOAD DATA INFILE '/docker-entrypoint-initdb.d/data/protein_info_clean.txt'
INTO TABLE protein_info
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(protein_id, protein_name, protein_size, annotation);
-- Add UniProt ID column
ALTER TABLE protein_info
ADD COLUMN uniprot_id VARCHAR(10) AFTER protein_name,
ADD COLUMN uniprot_name VARCHAR(16) AFTER uniprot_id;

-- DROP TABLE IF EXISTS tmp_uniprot_mapping;
CREATE TEMPORARY TABLE tmp_uniprot_mapping (
    uniprot_id VARCHAR(10),
    uniprot_name VARCHAR(16),
    protein_id VARCHAR(15),
    preferred_name VARCHAR(255)
);

-- LOAD DATA LOCAL INFILE "C:\\Users\\gotar\\GitHub\\BeltingLab\\Projects\\neo-surfaceome\\database\\uniprot_ID_mapping_clean.txt"
LOAD DATA INFILE '/docker-entrypoint-initdb.d/data/uniprot_ID_mapping_clean.txt'
INTO TABLE tmp_uniprot_mapping
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(uniprot_id, uniprot_name, protein_id, preferred_name);

UPDATE protein_info pi
JOIN tmp_uniprot_mapping um
  ON pi.protein_id = um.protein_id
SET pi.uniprot_id = um.uniprot_id,
    pi.uniprot_name = um.uniprot_name;
    
-- Protein info index
CREATE INDEX idx_protein_id ON protein_info(protein_id);

-- DROP TABLE IF EXISTS gene_info;
CREATE TABLE gene_info (
	transcript_id VARCHAR(15) PRIMARY KEY, 
    chromosome VARCHAR(5),
    _start INT,
    _stop INT,
    strand VARCHAR(1),
    _length INT,
    transcript_type VARCHAR(255),
    protein_id VARCHAR(15),
    gene_id VARCHAR(15),
    gene_symbol VARCHAR(15),
    hgnc_symbol VARCHAR(15)
);

-- LOAD DATA LOCAL INFILE "C:\\Users\\gotar\\GitHub\\BeltingLab\\Projects\\neo-surfaceome\\database\\gene_info.txt"
LOAD DATA INFILE '/docker-entrypoint-initdb.d/data/gene_info.txt'
INTO TABLE gene_info
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(transcript_id, chromosome, _start, _stop, strand, _length, transcript_type, protein_id, gene_id, gene_symbol, hgnc_symbol);

-- Gene info indexes
CREATE INDEX idx_transcript_id ON gene_info(transcript_id);
CREATE INDEX idx_gene_info_protein_id ON gene_info(protein_id);
CREATE INDEX idx_gene_info_gene_id ON gene_info(gene_id);

-- DROP TABLE IF EXISTS protein_seq;
CREATE TABLE protein_seq (
    protein_id VARCHAR(15),
    _length INT,
    sequence LONGTEXT
);

DROP TABLE IF EXISTS tmp_protein_seq;
CREATE TEMPORARY TABLE tmp_protein_seq (
    transcript_id VARCHAR(15),
    _length INT,
    sequence LONGTEXT
);

-- LOAD DATA LOCAL INFILE "C:\\Users\\gotar\\GitHub\\BeltingLab\\Projects\\neo-surfaceome\\database\\protein_sequences.txt"
LOAD DATA INFILE '/docker-entrypoint-initdb.d/data/protein_sequences.txt'
INTO TABLE tmp_protein_seq
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(transcript_id, _length, sequence);
-- Update transcript-based IDs to correct protein_ids
INSERT INTO protein_seq (protein_id, _length, sequence)
SELECT
    gi.protein_id,
    tps._length,
    tps.sequence
FROM tmp_protein_seq tps
JOIN gene_info gi ON tps.transcript_id = gi.transcript_id
WHERE gi.protein_id IS NOT NULL;
DROP TABLE IF EXISTS tmp_protein_seq;

-- DROP TABLE IF EXISTS transcript_seq;
CREATE TABLE transcript_seq (
	transcript_id VARCHAR(15) PRIMARY KEY, 
    _length INT,
    sequence LONGTEXT,
    FOREIGN KEY (transcript_id) REFERENCES gene_info(transcript_id)
);

-- LOAD DATA LOCAL INFILE "C:\\Users\\gotar\\GitHub\\BeltingLab\\Projects\\neo-surfaceome\\database\\transcript_sequences.txt"
LOAD DATA INFILE '/docker-entrypoint-initdb.d/data/transcript_sequences.txt'
INTO TABLE transcript_seq
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(transcript_id, _length, sequence);

-- DROP TABLE IF EXISTS protein_interactions;
CREATE TABLE protein_interactions (
	protein_A VARCHAR(15), 
    protein_B VARCHAR(15),
    neighborhood_score INT,
    fusion_score INT,
    cooccurence_score INT,
    coexpression_score INT,
    experimental_score INT,
    database_score INT,
    textmining_score INT,
    combined_score INT,
    FOREIGN KEY (protein_A) REFERENCES protein_info(protein_id),
    FOREIGN KEY (protein_B) REFERENCES protein_info(protein_id)
);

-- LOAD DATA LOCAL INFILE "C:\\Users\\gotar\\GitHub\\BeltingLab\\Projects\\neo-surfaceome\\database\\protein_links_clean.txt"
LOAD DATA INFILE '/docker-entrypoint-initdb.d/data/protein_links_clean.txt'
INTO TABLE protein_interactions
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(protein_A, protein_B, neighborhood_score, fusion_score, cooccurence_score,
 coexpression_score, experimental_score, database_score, textmining_score, combined_score
);
-- Protein Interactions
CREATE INDEX idx_protein_A ON protein_interactions(protein_A);
CREATE INDEX idx_protein_B ON protein_interactions(protein_B);

-- DROP TABLE IF EXISTS domain_info;
CREATE TABLE domain_info (
	domain_id VARCHAR(7) PRIMARY KEY, 
    domain_name VARCHAR(255),
    annotation VARCHAR(255)
);

-- LOAD DATA LOCAL INFILE "C:\\Users\\gotar\\GitHub\\BeltingLab\\Projects\\neo-surfaceome\\database\\pfam_domains.txt"
LOAD DATA INFILE '/docker-entrypoint-initdb.d/data/pfam_domains.txt'
INTO TABLE domain_info
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(domain_id, domain_name, annotation);
-- Domain info index
CREATE INDEX idx_domain_id ON domain_info(domain_id);

-- DROP TABLE IF EXISTS domain_annotation;
CREATE TABLE domain_annotation (
	transcript_id VARCHAR(15),
	protein_length INT,
    domain_id VARCHAR(7),
	domain_length INT,
	domain_start INT,
	domain_end INT,
    accuracy FLOAT,
	score FLOAT,
	E_value FLOAT,
    FOREIGN KEY (transcript_id) REFERENCES gene_info(transcript_id),
    FOREIGN KEY (domain_id) REFERENCES domain_info(domain_id)
);

LOAD DATA LOCAL INFILE "C:\\Users\\gotar\\GitHub\\BeltingLab\\Projects\\neo-surfaceome\\database\\transcript_domain_annotation.txt"
INTO TABLE domain_annotation
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(transcript_id, protein_length, domain_id, domain_length, domain_start, domain_end, accuracy, score, E_value);
-- Domain annotation indexes
CREATE INDEX idx_da_transcript_id ON domain_annotation(transcript_id);
CREATE INDEX idx_da_domain_id ON domain_annotation(domain_id);

-- DROP TABLE IF EXISTS domain_interactions;
CREATE TABLE domain_interactions (
	domain_A VARCHAR(15), 
    domain_B VARCHAR(15),
    FOREIGN KEY (domain_A) REFERENCES domain_info(domain_id),
    FOREIGN KEY (domain_B) REFERENCES domain_info(domain_id)
);

LOAD DATA LOCAL INFILE "C:\\Users\\gotar\\GitHub\\BeltingLab\\Projects\\neo-surfaceome\\database\\domain_interactions.txt"
INTO TABLE domain_interactions
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
(domain_A, domain_B);
-- Domain interactions indexes
CREATE INDEX idx_di_domain_A ON domain_interactions(domain_A);
CREATE INDEX idx_di_domain_B ON domain_interactions(domain_B);

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

CREATE TEMPORARY TABLE tmp_high_conf_protein_interactions AS
SELECT * 
FROM protein_interactions
WHERE combined_score >= 900;

CREATE TABLE protein_domain_map AS
SELECT gi.gene_id, pi.protein_id, di.domain_id, di.domain_name, da.domain_start, da.domain_end
FROM gene_info gi
JOIN protein_info pi ON gi.protein_id = pi.protein_id
JOIN tmp_high_conf_domain_annotation da ON gi.transcript_id = da.transcript_id
JOIN domain_info di ON da.domain_id = di.domain_id;

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