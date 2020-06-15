## USPTO to Neo4j

This package is meant to help data analysts and chemist who wish to view the USPTO data in Neo4j. Viewing US Patent
chemical data in Neo4j allows for chemical reaction pathway searching that traditional databases can not provide. 

## Getting required data

Make sure you have ~45 GB of free space before beginning this process. 

Due to the size of the USPTO data, it is not reasonable to package the data alongside of this package. Some additonal
steps are needed to download and prepare the data.

The raw data can be downloaded from
https://figshare.com/articles/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873.
It is important to download all the data, and as of writing this it is ~1.39 GB. It may be useful to save the zipped 
directory as a backup. Then uncompress everything, including the directories inside of the main folder. A prorgram such 
as 7-zip may be helpful in extracting the folders https://www.7-zip.org/download.html. Depending on operating systems, 
there may be different ways the folders get extracted. Make sure you end up with a file tree inside
of the directory similar to

* [5104873]()
    * [applications]()
        * [2001]()
        * [2002]()
        * [...]()
    * [grants]()
        * [1976]()
        * [1977]()
        * [...]()
        
Make sure the data is extracted and formatted properly, otherwise the initialization of the data will fail. 
When you download the data this way, the data is formatted only in a xml format. When running the initialization step in 
the next section, make sure to set ```convert_xml_to_csv=True```.

A more simplified way of downloading the data is going to the release for this package and downloading 5104873.7z from
https://github.com/andreshyer/USPTOtoNeo4j/releases/tag/v0.1-beta.
and simply download and extract the data. If you download the data this way, you must set
```convert_xml_to_csv=False```. 

## Inserting USPTO data into Neo4j

If you are completely unfamiliar with Neo4j, it is advised you learn some basics about Neo4j 
https://neo4j.com/docs/getting-started/current/.

The steps needed to insert data into Neo4j are:

1. Ensure python environment is confiured properly using conda and utn.yml
2. Create and startup a Neo4j graph database
3. Configure and run ```bulk_insert.py```

To configure ```bulk_insert.py```, it is necessary to hard code the parameters. The setting are found at the bottom of
```bulk_insert.py```. The defaults settings are
```
 if __name__ == "__main__":
     params = dict(
         patents_directory='/home/user/Desktop/5104873',
         number_of_cups=3,
         convert_xml_to_csv=False,
         clean_checker_files=False,
         insert_compounds_with_functional_groups=False,
         insert_change_in_functional_groups=False,
         loading_bars=True,
     )
 
     graph = Graph(user='neo4j', password='password', bolt_port='bolt://localhost:7687')
     InitPatentsIntoNeo4j(**params)     
```

Make sure the patents_directory refers the head of the directory tree as shown above.
The most important setting to note is ```clean_checker_files=False```. This clean checker files option keeps track
of which files have already been inserted into Neo4j. So when clean checker files is set to false, it is safe to pause
the script and resume later. But, if you want to create a new graph or re-initalize the graph, then set
```clean_checker_files=True```. This will delete all the tracking files and allow you to restart the inserting of the
files.

The other options are hopefully self explanatory, with the exception of ```insert_change_in_functional_groups=False```.
This option seeks to track what functional groups change between reactants and products in a given reaction. This option
is useful depending on certain tasks, but is set to false by default. 

Also, it is important to note that calculating functional groups will require much more time to insert the data into 
Neo4j. 

## Example queries

The schema is a very simple, but powerful, schema. There are two entities, compounds and reactions. And the
relationships are as follows:
* Compound-[reacts]->Reaction
* Compound-[catalyzes]->Reaction 
* Compound-[solvates]->Reaction 
* Reaction-[produces]->Compound
                
The following queries will use rdkit-canonical smiles. Smiles are perferred over chemcial names, as there are many
different chemical names for a compound, but only one canonical smiles string that represents it. So to ensure a 
match is made, smiles are used over names. Although it is possible to query using chemical names and inchi keys.    

Also, for these queries, it is assumed you have some working knowledge of the cypher language for Neo4j.             
                
So a simple query is to look for a reaction pathway between two compound, with the question 'Starting from (A), what
reactions are needed to get to (B)?'. The example query would become

```
MATCH path = (a:Compound {smiles: "Cc1cc([N+](=O)[O-])ccc1N"})-[:reacts|produces*..5]->(b:Compound {smiles: "Cc1cc([N+](=O)[O-])ccc1-n1c(=O)n(-c2ccc([N+](=O)[O-])cc2C)c(=O)n(-c2ccc([N+](=O)[O-])cc2C)c1=O"})
RETURN path LIMIT 5
```

Another query to ask is "Starting from (A), what compounds can be made that have a molar mass over 160 g/mol?"

```
MATCH path = (a:Compound {smiles: "O=Cc1cc(Br)cc(C=O)c1O"})-[:reacts|produces*..7]->(b:Compound)
WHERE b.smiles <> a.smiles`
AND b.molarmass > 160
RETURN path LIMIT 5
```                  

These queries lead to the very natural question "What compounds with a molar mass over 160 g/mol can 
be used to synthesize (B)?"

```
MATCH path = (a:Compound {smiles: "CC(C)(C)OC(=O)N1CCCC(c2nc(COc3ccc(-n4cnnn4)cc3)cs2)C1"})<-[:reacts|produces*..7]-(b:Compound)
WHERE b.smiles <> a.smiles
AND b.molarmass > 160
RETURN path LIMIT 5      
```     

The power of such a query is that if price data for compounds can ever be found. A simple replace of molar mass for 
price can allow the question "What is the cheapest starting material to synthesize (B)?" But, without price data, this
query is not possible. 