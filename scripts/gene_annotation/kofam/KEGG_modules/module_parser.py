#! /usr/bin/env python
"""
KEGG module parser compares KOs of a given module and KOs of a given species to 
determine if module is functional in given sample.

Usage:
python module_parser.py -1 kofam_output.tsv -2 list_of_kos_in_module.txt

Example: 
python module_parser.py -1 ./Plac_prot_ko_annotations.tsv -2 ./hist_biosyn.txt
"""


import csv
import re
import argparse

mp = argparse.ArgumentParser(description="Compares KOs of a given module and KOs of a given species to determine if module is functional in given sample")

mp.add_argument("-1", "--tsv_file", type=str, help="path to kofam output tsv file for desired species", required=True)
mp.add_argument("-2", "--module_file", type=str, help="path to txt file with module KOs", required=True)


def main():
    """
    Check if the KEGG module is functional given present KOs.
    """
    m = mp.parse_args()
    #Load files
    module_str = load_module_kos(m.module_file)
    present_kos = load_species_kos(m.tsv_file)

    #Retrieve boolean module KOs 
    module_KOs = parse_kegg_module(module_str)

    #Extract list of KOs present in module without separators 
    all_kos = extract_KOs(module_KOs)

    #Create boolean list of KOs
    safe_module_KOs = boolean_kos(present_kos, module_KOs)

    #Evaluates module boolean and identifying missing KOs
    functionality = eval_module(safe_module_KOs)

    #Identify missing KOs
    missing = [ko for ko in all_kos if ko not in present_kos]
    present = [ko for ko in all_kos if ko in present_kos]
    missing = " ".join(missing)
    present = " ".join(present)


    print(functionality, "\t", missing, "\t", present)



def load_species_kos(tsv_file: str) -> set[str]:
    """
    Reads a kofam output TSV file and extracts all KOs from the 3rd column into a set.
    """
    present_kos = set()
    with open(tsv_file, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        
        #Skip the headers
        next(reader, None)
        next(reader, None)
        
        for row in reader:
            if len(row) >= 3:
                ko = row[2].strip()
                if ko and ko.startswith('K') and ko[1:].isdigit():
                    present_kos.add(ko)
    return present_kos

                    

def load_module_kos(module_file: str) -> str:
    """
    Reads a text file including all the KOs for a specific module
    """
    with open(module_file, "r") as f:
        content = f.read()
        module_str = content.replace("\n", "")
    return module_str


#module_str = "K00801 K00511 K01852 K05917 (K00222,K19532) K07750 K07748 K13373 K09828 K01824 K00227 K00213"


#KOs = ["K00801","K00511", "K01852", "K05917", "(K00222,K19532)", etc]
#if "-k" in parenttheses then it will not see it ?

def parse_kegg_module(module_str: str) -> str:
    """
    Convert KEGG module notation into boolean expression
    """
    module_KOs = module_str.strip()

    # Split to skip optional KOs - doesnt work for -ko inside parethneses
    # KOs = module_KOs.split()
    # required_KOs = []
    # for KO in KOs:
    #     if KO.startswith('-K'):
    #         continue
    #     required_KOs.append(KO)
    # module_KOs = " ".join(required_KOs)

    #remove optional KOs
    module_KOs = re.sub(r'-K\d{5}', '', module_KOs)

    # "--" is KEGG's placeholder for a step with no assigned KO (a spontaneous
    # or uncharacterized reaction) — treat it as automatically satisfied by
    # dropping it entirely, rather than leaving it as an invalid operand in
    # the eval'd boolean expression. Matched only as a standalone token (not
    # preceded/followed by a non-space char) so it doesn't touch anything else.
    module_KOs = re.sub(r'(?<!\S)--(?!\S)', '', module_KOs)

    # Clean up empty parentheses like ()
    module_KOs = re.sub(r'\(\s*\)', '', module_KOs)
    # Normalize whitespace
    module_KOs = re.sub(r'\s+', ' ', module_KOs).strip()

    # Replace KEGG symbols with logical operators
    module_KOs = module_KOs.replace(' ', ' and ')
    module_KOs = module_KOs.replace(',', ' or ')
    module_KOs = module_KOs.replace('+', ' and ')

    return module_KOs

#"K00801 and K00511 and K01852 and

#module_KOs = "K00801 and K00511 and K01852 and  K05917 and (K00222 or K19532) and K07750 and K07748 and K13373 and K09828 and K01824 and K00227 and K00213"
#K00232 and K12405 and (K07513 or K08764)

def extract_KOs(module_KOs:str) -> list:
    """
    Extract all KOs from the module and store as a list 
    """
    all_kos = re.findall(r'K\d{5}', module_KOs)
    return all_kos

    # Split on characters that separate KOs and replace with spaces to split
    # separators = [' and ', ' or ', ' and (', ') and ', ' or (', ') or ']
    # for sep in separators:
    #     module_KOs = module_KOs.replace(sep, ' ')
    # module_KOs = module_KOs.split()
    

    # Keep only valid KO IDs (K followed by 5 digits) and remove duplicates with set
    # kos = []
    # for ko in module_KOs:
    #     if len(ko) == 6 and ko.startswith('K') and ko[1:].isdigit():
    #         kos.append(ko)
    # all_kos = list(set(kos))
    

#['K00232', 'K12405']

def boolean_kos(present_kos:str, module_KOs:str) -> str:
    """
    Replace each KO with True or False depending on whether it's in present_kos
    """
    found_kos = re.findall(r'K\d{5}', module_KOs) #get around parentheses problem
    for ko in found_kos:
        if ko in present_kos:
            module_KOs = module_KOs.replace(ko, "True")
        else:
            module_KOs = module_KOs.replace(ko, "False")

    return module_KOs
    

#safe_modlue_KOs = "True and True and (True or False) and True"



def eval_module(safe_module_KOs: str,) -> bool:
   """ 
   Evaluates module boolean and identifying missing KOs.
   If every required KO (and every required group) is satisfied then functional = True
   """
   #eval
   try:
    functionality = eval(safe_module_KOs)
   except Exception as e:
    print("Error evaluating:", safe_module_KOs)
    raise e
   return functionality




if __name__ == "__main__":
    main()



# (True) and (True and True or True or True) and (True or True) and (True and True or True)((True and True and (False or True or True)) or (True and True)) and (True or True)
# (False) and (False and False or False or True) and (True or True) and (True and True or True)((False and True and (False or True or True)) or (True and True)) and (False or True)
# (K00765-K02502) (K01523 K01496,K11755,K14152) (K01814,K24017) (K02501+K02500,K01663) ((K01693 K00817 (K04486,K05602,K18649)),(K01089 K00817)) (K00013,K14152)