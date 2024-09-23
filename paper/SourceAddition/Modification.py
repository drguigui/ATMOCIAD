#!/usr/bin/env python3
import re

# Function to read bibtex file and extract citation keys
def extract_citation_keys(bibtex_file):
    with open(bibtex_file, 'r') as file:
        bibtex_content = file.read()
    
    # Regular expression to capture citation keys like Tobiska1999
    citation_keys = re.findall(r'@\w+\{([^,]+),', bibtex_content)
    return citation_keys

# Function to replace "Author Year" with "\cite{AuthorYear}" in the tex document
def replace_author_year_with_citations(tex_file, output_file, citation_keys):
    with open(tex_file, 'r') as file:
        tex_content = file.read()
    
    # Build a case-insensitive regex for each citation key
    for key in citation_keys:
        # Extract the author name and year from the citation key (e.g., Tobiska1999 -> Tobiska 1999)
        match = re.match(r'([A-Za-z]+)(\d{4})', key)
        if match:
            author, year = match.groups()
            
            # Create a regex pattern for "Author Year" case-insensitive match
            pattern = re.compile(rf' {author}\s+{year}', re.IGNORECASE)
            
            # Replace occurrences of "Author Year" with "\cite{AuthorYear}"
            tex_content = pattern.sub(' \\\cite{'+key+'}', tex_content)

        match = re.match(r'([A-Za-z]+)(\d{4})', key)
        if match:
            author, year = match.groups()
          
           # Create a regex pattern for "Author et al. Year" case-insensitive match and replace it with "\cite{AuthorYear}"
            pattern = re.compile(rf' {author}\s+et\s+al.\s+{year}', re.IGNORECASE)
            tex_content = pattern.sub(' \\\cite{'+key+'}', tex_content)
            pattern = re.compile(rf' {author}\s+et\s+al\s+{year}', re.IGNORECASE)
#           pattern = re.compile(rf'{author}\s+{year}', re.IGNORECASE)
           
           # Replace occurrences of "Author Year" with "\cite{AuthorYear}"
            tex_content = pattern.sub(' \\\cite{'+key+'}', tex_content)


    # We want to replace Avakyan 98 by \cite{Avakyan1998}
    pattern = re.compile(rf' Avakyan \d\d', re.IGNORECASE)
    tex_content = pattern.sub(r'\\cite{Avakyan1998}', tex_content)
    # We want to replace Avakyan by \cite{Avakyan1998}
    pattern = re.compile(rf' Avakyan ', re.IGNORECASE)
    tex_content = pattern.sub(r'\\cite{Avakyan1998}', tex_content)
    pattern = re.compile(rf'Gentieu and Mentall \d\d\d\d ', re.IGNORECASE)
    tex_content = pattern.sub( r'\\cite{Gentieu1973}', tex_content)
    pattern = re.compile(rf'Glass-Maujean and Schmoranzer \d\d\d\d ', re.IGNORECASE)
    tex_content = pattern.sub( r'\\cite{Glass2005}', tex_content)
    pattern = re.compile(rf' Singhal ', re.IGNORECASE)
    tex_content = pattern.sub(r'\\cite{Singhal2009}', tex_content)


    # Save the modified content to a new file
    with open(output_file, 'w') as file:
        file.write(tex_content)
def Change(bibtex_file, tex_directory):
    # We save the current directory
    current_directory = os.getcwd()


    # Extract citation keys from the BibTeX file
    citation_keys = extract_citation_keys(bibtex_file)
    os.chdir(tex_directory)
    save_directory = "./save"
    os.makedirs(save_directory, exist_ok=True)

    # We are going through all the tex files in the directory
    for file in glob.glob("*.tex"):
        tex_file = file
        output_file = os.path.join(save_directory, file)
        # Replace "Author Year" in the tex file with citations
        replace_author_year_with_citations(tex_file, output_file, citation_keys)
        print("Citations replaced successfully!")

    os.system("mkdir save2")
    for file in glob.glob("*.tex"):
        os.system("cp "+ file+" save2")
    # Now we are going to swap the files 
    os.chdir(save_directory)
    for file in glob.glob("*.tex"):
        os.system("cp "+ file+" ../")
    os.chdir(current_directory)
    
    print("Citations replaced successfully!")

# Main logic to extract citations and replace them in tex documents

import glob
import os
if __name__ == "__main__":
    tex_directory = "../../electron/resultat/"
    bibtex_file = '../atmociad.bib'  # Path to your bibtex file
    Change(bibtex_file, tex_directory)
    tex_directory = "../../photon/resultat/"
    Change(bibtex_file, tex_directory)
    tex_directory = "../../hydrogen/resultat/"
    Change(bibtex_file, tex_directory)
    tex_directory = "../../proton/resultat/"
    Change(bibtex_file, tex_directory)



