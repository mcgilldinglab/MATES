To download gene and TE (Transposable Element) references for your preferred species from the UCSC Genome Browser, follow these steps:

1. Go to the UCSC Genome Browser website: [https://genome.ucsc.edu/cgi-bin/hgGateway](https://genome.ucsc.edu/cgi-bin/hgGateway).

2. In the search bar, type the name of your preferred species and select it from the search results.

<img width="1512" alt="Screenshot 2024-11-06 at 16 56 31" src="https://github.com/user-attachments/assets/883a8948-6501-4409-894a-fc96f38e1f50">


3. Once you've selected your species, click on **Tools** in the upper tab.

4. From the dropdown menu, choose **Table Browser**.

<img width="1512" alt="Screenshot 2024-11-06 at 16 56 57" src="https://github.com/user-attachments/assets/4935ccc2-f662-4a7d-bf9a-4b66c22f0bb5">



5. You will be directed to the Table Browser. This step initializes the **clade** and **genome** fields for your selected species.

6. Now, follow these steps for Gene Reference:

   - In the **clade** dropdown menu, ensure the clade is set to the appropriate classification.
   - In the **genome** dropdown menu, select the genome assembly you want to use.
   - In the **Group** dropdown menu, select the one related to genes.
   - Set the **output format** as GTF.
   - Click the **'bigZip/genes'** to get formatted gene identifiers. 
   - The remaining fields can be set to their default values unless you have specific requirements.

<img width="1512" alt="Screenshot 2024-11-06 at 16 57 34" src="https://github.com/user-attachments/assets/1113efe7-e274-4a94-a1b0-f8ed6bc48bb1">

After opening the webpage, you can use either *ensGene.gtf* or *ncbiRefSeq.gtf* file.

<img width="1512" alt="Screenshot 2024-11-06 at 16 57 45" src="https://github.com/user-attachments/assets/6fb8b067-947b-4892-8874-5f3b11790796">

## Alternatively, you can use [Ensembl](https://useast.ensembl.org/) or related websites like [EnsemblPlants](https://plants.ensembl.org/index.html).

NCBI database does not always have the GTF with gene identifiers, such as Arabdopsis. We can use the Ensembl/EnsemblPlants/EnsemblFungi/../.

### Arabdopsis in EnsemblPlants

Download the GFF3 file.

<img width="1512" alt="Screenshot 2024-11-06 at 15 57 49" src="https://github.com/user-attachments/assets/7fc6df0f-ca91-40b2-a90b-69ee90b4be96">

### Drosophila

<img width="1511" alt="Screenshot 2024-11-06 at 16 00 17" src="https://github.com/user-attachments/assets/09aa4b43-33bc-4148-8046-46167457fde2">

You can use either GTF or GFF3 file.

<img width="1512" alt="Screenshot 2024-11-06 at 16 00 31" src="https://github.com/user-attachments/assets/a9d995d4-b206-43c8-a8e1-69509d4cc416">


7. For TE (Transposable Element) Reference:

   - Follow the same procedure as for the Gene Reference to initialize the **clade** and **genome** fields.
   - In the **Group** dropdown, you can use either *'all tracks'* or *repeat* related one.
   - In the **track** dropdown, select *RepeatMasker*.
   - Change the remaining fields as shown below, making sure they match what you used for the Gene Reference unless you have special requirements.

<img width="1506" alt="Screenshot 2024-11-06 at 17 18 47" src="https://github.com/user-attachments/assets/19c2e631-b295-4108-8240-5e27f15611c3">


   - Click the **get output** button.

8. After clicking the **get output** button for either the Gene Reference or TE Reference, the respective reference data will be processed, and you will be prompted to download the zipped reference file.

Follow these steps for both the Gene and TE References to obtain the required reference files for your preferred species.

Here is the example script to build reference genome after downloading the data.

```bash
gzip -d dm6.ensGene.gtf.gz
python build_reference.py --species Other --other_species_TE Drosophila_TE.csv --other_species_GTF dm6.ensGene.gtf
```
