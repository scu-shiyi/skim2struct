import re
from ete3 import Tree
class TreeLoad:
    def __init__(self, tree_file):

        self.file = tree_file
        self.condensed_tree = "" 
        self.newick_tree = ""
        self.species_data = []  
        self.raw_content = []  
        self.translate_dict = {}
        self.read_tree_file()


    def read_tree_file(self):

        try:
            with open(self.file, "r", encoding="utf-8") as file:
                self.raw_content = [line.strip() for line in file.readlines()]
        except Exception as e:
            print(f"Failed to read tree file: {e}")

    def parse_tree(self):

        species_map = {}
        species_id_counter = 1
        in_translate_block = False
        has_translate = False
        for line in self.raw_content:

            if not line:
                continue
   
            if "TRANSLATE" in line.upper():
                has_translate = True 
                in_translate_block = True
                continue
            if in_translate_block:
                if ";" in line:
                    in_translate_block = False
                    continue
                if line:
                    parts = line.replace(";", "").split()
                    if len(parts) >= 2:
                        species_id, species_name = parts[0], parts[1].replace("'", "").strip(",")
                        self.translate_dict[species_id] = species_name
                        # self.reverse_translate_dict['speices_name'] = species_name
                        self.species_data.append((species_id, species_name))
                    continue

            if "(" in line and ")" in line and ":" in line and ";" in line:
                line = re.sub(r"\[&&NHX:[^\]]+\]", "", line)  
                line = re.sub(r"\[&[^\]]*\]", "", line)
                line = re.sub(r"^\s*tree\s+\w+\s*=\s*", "", line, flags=re.IGNORECASE)  
                current_tree = line.strip()
                if has_translate: # nexus

                    for id, name in self.translate_dict.items():
                        current_tree = re.sub(rf'\b{id}', name, current_tree)
                    self.newick_tree = current_tree
                else: 
                    self.newick_tree = current_tree
                    matches = re.findall(r"([A-Za-z0-9][A-Za-z0-9._-]*):[\d.]+", current_tree)
                    for species_name in matches:
                        if species_name not in species_map:
                            species_map[species_name] = str(species_id_counter)
                            self.species_data.append((str(species_id_counter), species_name))
                            species_id_counter+=1

            species =set([i[1] for i in self.species_data])

        
        return species, self.newick_tree



