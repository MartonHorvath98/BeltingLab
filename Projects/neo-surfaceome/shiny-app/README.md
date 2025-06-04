# NSFW, Neo-Surfaceome Feature Workbench

- [NSFW, Neo-Surfaceome Feature Workbench](#nsfw-neo-surfaceome-feature-workbench)
  - [1. Application overview](#1-application-overview)
  - [2. Application configurations](#2-application-configurations)
  - [3. Functionality](#3-functionality)
  - [4. Dockerization](#4-dockerization)

  This [R Shiny](https://shiny.posit.co/) application provides an interactive graphical user interface (GUI) for querying and visualizing variant-related domain disruptions in surface proteins. It connects to our [MySQL database](../database/README.md) high-confidence, isoform-resolved protein–protein interaction (PPI) data to support precision oncology insights in the future... :dizzy: :rocket:
  
  > [!NOTE]
  > *...for now, it only exists as a proof-of-concept execution to demonstrate the working connection between the components and the queriebility of the database.*

## 1. Application overview

### Main features:
- Searchable table of the SURFME catalogue, surfaceome annotated genes and proteins
- Domain annotation viewer for selected genes
- Perturbed domain predictions from variant input (SNP position)

## 2. Application configuration

The application runs on localhost:
```yaml
host: 0.0.0.0
port: 8180
```
Which is explicitely set in the `app.R` file:
```r
options(shiny.host = "0.0.0.0")
options(shiny.port = 8180)
```
The `RMySQL` package is used to set up connection with the MySQL database, to forward user queries and return resulting tables. When the app server starts, a connection is set up with the MySQL database:
```r
SFcon <- dbConnect(
    RMySQL::MySQL(),
    dbname = db_name,
    host = db_host,
    port = db_port,
    user = db_user,
    password = db_password
  )
```
Where the following environment variables configure the connection:

| Variable | Default value | Description |
| --- | --- | --- |
| `DB_HOST` |  `127.0.0.1` | MySQL host (in Docker: `db`) |
| `DB_USER` |  `root` | MySQL user |
| `DB_PASSWORD` |  ***(must be defined by the user)*** | MySQL server password |
| `DB_NAME` |  `InteractomeDB` | MySQL database name |
| `DB_PORT` |  `root` | MySQL port |

## 3. Functionality

### Default SURFME gene table view

When the application is started, it connects to the MySQL server using the credential set via environment variable provided at runtime (e.g., in `compose.yml`). Then, it executes a JOIN query across `gene_info`, `surfme_filter`, and `protein_info` tables and merges results with the more detailed local annotation file of the SURFME proteins from `meta/SURFME_v2023.xlsx` on the uniprot_id column. On the GUI an interactive DT::DataTable is displayed with the columns: `UniprotID`, `GeneID`, `Symbol` and `Annotation`.

*e.g., first 5 rows of the table; print_data*
| | UniprotID | GeneID | Symbol | Annotation |
| :-: | :-: | :-: | :-: | --- |
| 1 | O43657 | ENSG00000000003 | TSPAN6 | Tetraspanin-6 (Tspan-6) (A15 homolog) (Putative NF-kappa-B-activating protein 321) (T245 protein) (Tetraspanin TM4-D) (Transmembrane 4 superfamily member 6) |
| 2 | Q6P499 | ENSG00000001461 | NIPAL3 | NIPA-like protein 3 |
| 3 | Q9Y6X5 | ENSG00000001561 | ENPP4 | Bis(5'-adenosyl)-triphosphatase ENPP4 (EC 3.6.1.29) (AP3A hydrolase) (AP3Aase) (Ectonucleotide pyrophosphatase/phosphodiesterase family member 4) (E-NPP 4) (NPP-4) |
| 4 | P13569 | ENSG00000001626 | CFTR | Cystic fibrosis transmembrane conductance regulator (CFTR) (ATP-binding cassette sub-family C member 7) (Channel conductance-controlling ATPase) (EC 5.6.1.6) (cAMP-dependent chloride channel) |
| 5 | P14209 | ENSG00000002586 | CD99 | CD99 antigen (12E7) (E2 antigen) (Protein MIC2) (T-cell surface glycoprotein E2) (CD antigen CD99) |

### Domain annotation view

Individual rows can be selected from the DataTable:
```r
DT::datatable(
    print_data, 
    selection = "single",
    options = list(scrollX = TRUE, scrollY = "300px", paging = FALSE, 
                   dom = '<"top"lf>rt<"bottom"><"clear">'),
    style = "bootstrap",
    class = "compact stripe hover row-border order-column"
)
```
When a row is selected, the app runs a secondary query on the `protein_domain_map` table via the selected GeneID and then displays unique, best-ranked domain annotations (lowest E-value) per domain ID. The resulting domain table is displayed as a non-reactive block, only showing the values in a clean table without pagination or filter options:
```r
domain_query <- paste0("
    SELECT DISTINCT pdm.domain_id, pdm.domain_name, pdm.domain_start, pdm.domain_end, pdm.accuracy, pdm.E_value
    FROM gene_info gi
    JOIN protein_domain_map pdm ON gi.gene_id = pdm.gene_id
    WHERE gi.gene_id = '", selected_gene_id, "';")
      
domain_data <- dbGetQuery(SFcon, domain_query, params = list(selected_gene_id))
domain_data <- domain_data %>% 
    dplyr::group_by(domain_id) %>% 
    dplyr::arrange(E_value, .by_group = TRUE) %>%
    dplyr::slice(1)
      
# Render as DT if any domains found
if (nrow(domain_data) > 0) {
    # Save domain_data to reactive value
    output$domainTable <- DT::renderDataTable({
        DT::datatable(
            domain_data,
            rownames = FALSE,
            options = list(
                paging = FALSE,      # no pagination
                searching = FALSE,   # no search box
                info = FALSE,        # no 'Showing 1 to n of n entries'
                ordering = FALSE,    # no sorting
                dom = 't'            # only the table
            ),
            style = "bootstrap",
            class = "compact stripe hover",
            escape = FALSE
        )
    })
}
```
> [!NOTE]
> *This was meant to be the stored procedure `MapDomains` although a little touch-up was needed. It will be corrected...*

### Detection of perturbed domain-domain interactions

Clicking the `“PerturbedDomain”` button reveals two inputs: a `textInput` (GeneID) and `numericInput` (SNP Position) fields, and an additional action button. For now - while the application is in the current proof-of-concept phase - the user can manually submit values to call the stored procedure with:
```sql
CALL PerturbedDomains('GENEID', POSITION);
```
The query outputs the affected domain(s) caused by the variant location, replacing the default SURFME gene table view the `Perturbed Domains` table, including four columns: `affected gene`, `affected doamin`, `potentially disrupted ddi` and `potentially disrupted partner`. *e.g., "CALL PeturbedDomains("ENSG00000146648", 100);"*
| affected_gene | affected_domain | potentially_disrupted_ddi | potentially_disrupted_partner |
| :-: | :-: | :-: | :-: |
| EGFR | PF07714 | PF07714 | PTK2 |
| EGFR | PF07714 | PF18038 | PTK2 |
| EGFR | PF07714 | PF00017 | RASA1 |
| EGFR | PF07714 | PF00018 | RASA1 |
| EGFR | PF07714 | PF00307 | VAV2 |
| EGFR | PF07714 | PF00018 | VAV2 |
| EGFR | PF07714 | PF14604 | VAV2 |
| EGFR | PF07714 | PF00017 | VAV2 |
| EGFR | PF07714 | PF00130 | VAV2 |
| EGFR | PF07714 | PF00169 | VAV2 |


> [!IMPORTANT]
> ***QC has not been realized yet on he user input. Non-existent 'geneID'-s or out-of-range 'SNV positions' do not break the program, but the reponse is not curated to explain the user's error.***

# 4. Dockerization

The app is containerized using rocker/shiny as base:
```dockerfile
FROM rocker/shiny
```
In the Docker build the following additional key steps are executed:
- installing system libraries for MySQL:
```bash
apt-get install -y default-libmysqlclient-dev
```
- installing R packages:
```r
install.packages(c('shiny', 'DBI', 'RMySQL', 'DT', 'dplyr', 'openxlsx'))
```
- setting up directory structure:
```bash
/home/shiny-app/
├── app.R
├── meta/SURFME_v2023.xlsx
└── www/nsfw-logo.png
```
- exposing port `8180` and runing the `shiny::runApp()` on container start.

### Docker Compose Setup
```yaml
version: "3.8"
# Docker Compose file for setting up a MySQL database and a Shiny app
services:
  db:
    image: mysql:8.0
    container_name: neo-db
    restart: always
    environment:
      MYSQL_ROOT_PASSWORD: pw # replace with actual pasword
      MYSQL_DATABASE: InteractomeDB
    volumes:
      - ./database/init.sql:/docker-entrypoint-initdb.d/init.sql
      - ./database/data:/docker-entrypoint-initdb.d/data
    command: --local-infile=1
    ports:
      - "3307:3306"
    networks:
      - neo-net

  shiny:
    build:
      context: ./shiny-app
    container_name: neo-shiny
    ports:
      - "8180:8180"
    depends_on:
      - db
    environment:
      DB_HOST: db
      DB_USER: root
      DB_PASSWORD: pw # replace with actual pasword
      DB_NAME: InteractomeDB
      DB_PORT: 3306
    networks:
      - neo-net

networks:
  neo-net:
```
- the database (neo-db) is started with init scripts and local data.
- the Shiny app (neo-shiny) connects to the db service via internal Docker networking.
- the app listens on localhost:8180 (accessible in browser).

### Folder structure
```swift
project/
├── shiny-app/
│   ├── app.R
│   ├── meta/
│   │   └── SURFME_v2023.xlsx
│   └── www/
│       └── nsfw-logo.png
├── database/
│   ├── init.sql
│   └── data/
│       └── ...  # .txt imports
├── docker-compose.yml
└── README.md
```

:rewind: *[Return](../README.md) to the main README file...* 
