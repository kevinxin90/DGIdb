import hashlib
import os
import requests

from biothings.utils.dataload import tabfile_feeder, dict_sweep, unlist
import yaml

"""
  | Column Names             | Key Names
--+--------------------------+---------------------------
0 | gene_name                | subject.SYMBOL
1 | gene_claim_name          |
2 | entrez_id                | subject.NCBIGene
3 | interaction_claim_source | association.provided_by
4 | interaction_types        | association.relation_name
5 | drug_claim_name          |
6 | drug_claim_primary_name  |
7 | drug_name                | object.name
8 | drug_chembl_id           | object.CHEMBL_COMPOUND
9 | PMIDs                    | association.pubmed
"""


def load_annotations(data_folder):
    data_file = os.path.join(data_folder, "interactions.tsv")
    data = tabfile_feeder(data_file, header=0)
    header = next(data)

    mapfile = open(os.path.join(
        data_folder, "predicate-remap.yaml"), 'r').read()
    remappingdata = yaml.safe_load(mapfile)

    def get_predicate(relation):
        if relation != "":
            key = ':'.join(("DGIdb", relation))
            return remappingdata[key]['rename'][0]
        return ""

    def get_gene_id(gene_name):
        query = ("http://mygene.info/v3/query?q=symbol:{}"
                 "&fields=entrezgene&species=human".format(gene_name))
        response = requests.get(query)
        if response.status_code == "200":
            data = response.json()
            entrez_id = data['hits'][0]['entrezgene']
            if entrez_id != "":
                return entrez_id
        return None

    def get_chem_id(drug_name):
        query = ("http://mychem.info/v1/query?q=chembl.pref_name:{}"
                 "&fields=chembl.molecule_chembl_id".format(drug_name))
        response = requests.get(query)
        if response.status_code == "200":
            data = response.json()
            chembl_id = data['hits'][0]['chembl']['molecule_chembl_id']
            if chembl_id != "":
                return chembl_id
        return None

    for rec in data:
        # Create a hash for _id
        bytestr = bytearray("-".join(rec), 'utf-8')
        hashstr = hashlib.blake2b(bytestr, digest_size=8).hexdigest()

        # Document framework
        doc = {
            "_id": hashstr,
            "subject": {},
            "object": {},
            "association": {}
        }

        # Subject
        entrez_id = rec[header.index("entrez_id")]
        gene_name = rec[header.index("gene_name")]
        if entrez_id == "":
            if gene_name == "":
                continue  # Skip the record
            resp = get_gene_id(gene_name)
            if resp is None:
                subject_id = 'name:' + gene_name
            else:
                entrez_id = resp
                subject_id = 'NCBIGene:' + resp
        else:
            subject_id = 'NCBIGene:' + entrez_id
        doc['subject']['NCBIGene'] = entrez_id
        doc['subject']['SYMBOL'] = gene_name
        doc['subject']['id'] = subject_id

        # Object
        drug_name = rec[header.index("drug_name")]
        drug_chembl_id = rec[header.index("drug_concept_id")]
        if drug_chembl_id == "":
            if drug_name == "":
                continue  # Skip the record
            resp = get_chem_id(drug_name)
            if resp is None:
                object_id = 'name:' + drug_name
            else:
                drug_chembl_id = resp
                object_id = 'CHEMBL.COMPOUND:' + resp
        elif drug_chembl_id.startswith("chembl:"):
            object_id = 'CHEMBL.COMPOUND:' + drug_chembl_id.split(':')[-1]
            drug_chembl_id = 'CHEMBL.COMPOUND:' + drug_chembl_id.split(':')[-1]
        doc['object']['name'] = drug_name
        doc['object']['CHEMBL_COMPOUND'] = drug_chembl_id
        doc['object']['id'] = object_id

        # Association
        interaction_types = rec[header.index("interaction_types")].replace(
            " ", "_").split(",")
        pmids = rec[header.index("PMIDs")].split(",")
        interaction_claim_source = rec[header.index(
            "interaction_claim_source")]
        edge_labels = []
        for interaction in interaction_types:
            edge_labels.append(get_predicate(interaction))
        if edge_labels == ['']:
            edge_labels = 'physically_interacts_with'
        doc['association']['edge_label'] = edge_labels
        doc['association']['relation_name'] = interaction_types
        doc['association']['pubmed'] = pmids
        doc['association']['provided_by'] = interaction_claim_source
        doc['association']['interaction_group_score'] = rec[header.index(
            "interaction_group_score")]
        # Cleanup
        doc = dict_sweep(doc)
        doc = unlist(doc)
        yield doc


if __name__ == "__main__":
    import json

    annotations = load_annotations("data/")
    for a in annotations:
        print(json.dumps(a, indent=2))
