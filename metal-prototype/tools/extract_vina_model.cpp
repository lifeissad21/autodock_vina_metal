#include <iomanip>
#include <iostream>
#include <string>

#include "parse_pdbqt.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: extract_vina_model ligand.pdbqt\n";
        return 2;
    }
    model parsed = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
    const ligand lig = parsed.get_ligand(0);
    std::cout << std::setprecision(9);
    for (std::size_t i = 0; i < parsed.num_movable_atoms(); ++i) {
        const atom atom = parsed.get_atom(i);
        const vec coordinates = parsed.get_coords(i);
        std::cout << "ATOM " << i << ' ' << atom.xs << ' '
                  << coordinates[0] << ' ' << coordinates[1] << ' ' << coordinates[2] << '\n';
    }
    for (const interacting_pair& pair : lig.pairs) {
        std::cout << "PAIR " << pair.a << ' ' << pair.b << '\n';
    }
    std::cout << "TORSIONS " << parsed.get_size().ligands.front() << '\n';
    return 0;
}
