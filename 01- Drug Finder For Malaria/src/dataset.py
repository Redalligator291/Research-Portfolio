import mysql.connector
from mysql.connector import Error

try:
    # Establish connection to the MySQL database
    connection = mysql.connector.connect(
        host='localhost',
        database='chembl_db',
        user='root',
        password='root'
    )

    if connection.is_connected():
        cursor = connection.cursor()

        # SQL query to extract drug-like compounds
        query = """
        SELECT 
            molecule_dictionary.chembl_id,
            compound_structures.canonical_smiles,
            compound_properties.full_mwt,
            compound_properties.alogp,
            compound_properties.hba,
            compound_properties.hbd,
            compound_properties.psa,
            compound_properties.rtb
        FROM 
            molecule_dictionary
        JOIN 
            compound_structures ON molecule_dictionary.molregno = compound_structures.molregno
        JOIN 
            compound_properties ON molecule_dictionary.molregno = compound_properties.molregno
        WHERE 
            compound_properties.full_mwt <= 500
            AND compound_properties.alogp <= 5
            AND compound_properties.hba <= 10
            AND compound_properties.hbd <= 5
            AND compound_properties.psa <= 140
            AND compound_properties.rtb <= 10
        """

        cursor.execute(query)
        results = cursor.fetchall()

        # Process and save the results
        with open('chembl_drug_like_compounds.csv', 'w') as f:
            f.write('chembl_id,smiles,molecular_weight,alogp,hba,hbd,psa,rtb\n')
            for row in results:
                f.write(','.join(map(str, row)) + '\n')

        print(f"Extracted {len(results)} compounds and saved to chembl_drug_like_compounds.csv")

except Error as e:
    print(f"Error connecting to MySQL database: {e}")

finally:
    if connection.is_connected():
        cursor.close()
        connection.close()
        print("MySQL connection is closed")