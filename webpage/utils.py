_GENETIC_CODES = [1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33]
_GENETIC_CODES_MITO = [2,3,4,5,9,13,14,16,21,22,23,24,33]

gencodes = {
    1: '1. The Standard Code',
    2: '2. The Vertebrate Mitochondrial Code',
    3: '3. The Yeast Mitochondrial Code',
    4: '4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code',
    5: '5. The Invertebrate Mitochondrial Code',
    6: '6. The Ciliate, Dasycladacean and Hexamita Nuclear Code',
    9: '9. The Echinoderm and Flatworm Mitochondrial Code',
    10: '10. The Euplotid Nuclear Code',
    11: '11. The Bacterial, Archaeal and Plant Plastid Code',
    12: '12. The Alternative Yeast Nuclear Code',
    13: '13. The Ascidian Mitochondrial Code',
    14: '14. The Alternative Flatworm Mitochondrial Code',
    15: '15. Blepharisma Nuclear Code',
    16: '16. Chlorophycean Mitochondrial Code',
    21: '21. Trematode Mitochondrial Code',
    22: '22. Scenedesmus obliquus Mitochondrial Code',
    23: '23. Thraustochytrium Mitochondrial Code',
    24: '24. Rhabdopleuridae Mitochondrial Code',
    25: '25. Candidate Division SR1 and Gracilibacteria Code',
    26: '26. Pachysolen tannophilus Nuclear Code',
    27: '27. Karyorelict Nuclear Code',
    28: '28. Condylostoma Nuclear Code',
    29: '29. Mesodinium Nuclear Code',
    30: '30. Peritrich Nuclear Code',
    31: '31. Blastocrithidia Nuclear Code',
    33: '33. Cephalodiscidae Mitochondrial UAA-Tyr Code',
}
GENETIC_CODES = list(gencodes.keys())
GENETIC_CODES_MITO = [k for k,v in gencodes.items() if "Mito" in v]

assert GENETIC_CODES == _GENETIC_CODES
assert GENETIC_CODES_MITO == _GENETIC_CODES_MITO


def gencode_id2title(s: str):
    return gencodes.get(s, None)

