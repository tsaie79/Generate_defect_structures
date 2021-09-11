__author__ = "jengyuantsai"

from pycdt.core.defectsmaker import ChargedDefectsStructures

from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.defects.core import Vacancy, Substitution
from pymatgen.io.vasp.inputs import Structure

import os
import numpy as np

class Tools:
    @classmethod
    def set_vacuum(cls, orig_st, vacuum):
        from mpinterfaces.utils import ensure_vacuum
        st = ensure_vacuum(orig_st, vacuum)
        return st

    @classmethod
    def phonopy_structure(cls, orig_st):
        """
        Phonopy and pos2aBR (https://github.com/zjwang11/UnconvMat) are to be installed
        :param orig_st: host material
        :return: standardize host material, pos2aBR outputs
        """
        from subprocess import call, check_output
        import shutil
        os.makedirs("standardize_st", exist_ok=True)
        os.chdir("standardize_st")
        orig_st.to("poscar", "POSCAR")
        call("phonopy --symmetry --tolerance 0.01 -c POSCAR".split(" "))
        std_st = Structure.from_file("PPOSCAR")
        std_st.to("poscar", "POSCAR")
        pos2aBR_out = check_output(["pos2aBR"], universal_newlines=True).split("\n")
        std_st = Structure.from_file("POSCAR_std")
        os.chdir("..")
        shutil.rmtree("standardize_st")
        return std_st, pos2aBR_out

    @classmethod
    def get_rand_vec(cls, distance): #was 0.001
        # deals with zero vectors.
        vector = np.random.randn(3)
        vnorm = np.linalg.norm(vector)
        return vector / vnorm * distance if vnorm != 0 else cls.get_rand_vec(distance)

    @classmethod
    def find_scaling_for_2d_defect(cls, pc, min_lc=15):
        x = 1
        while True:
            st = pc.copy()
            st.make_supercell([x,1,1])
            if st.lattice.a >= min_lc:
                break
            else:
                x += 1

        y = 1
        while True:
            st = pc.copy()
            st.make_supercell([1, y, 1])
            if st.lattice.b >= min_lc:
                break
            else:
                y += 1

        scaling = [x, y, 1]
        st = pc.copy()
        st.make_supercell(scaling)
        print("scaling matrix: {}".format(scaling))
        return scaling, st




class GenDefect:
    """
    Generate defect structures including vacancies and substitutions.
    Combining with some features including adding vacuum, distortion, and standraizing structure.
    """
    def __init__(self, orig_st, defect_type, natom, vacuum_thickness=None, distort=None, sub_on_side=None, standardize_st=True):
        """

        :param orig_st: host material (pymatgen Structure)
        :param defect_type:(defect_type(str), defect_index(int))
            defect = ChargedDefectsStructures(pc, antisites_flag=True).defects; defect[defect_type][defect_index]["name"]
        :param natom: number of atoms of supercell without defect
        :param vacuum_thickness: [vacuum thickness (float)] (list of float)
        :param distort: randomly move defct and its adjacent atoms by distance Å (float)
        :param sub_on_side: make substitution of one of the adjacent atoms (list of str e.g. ["Si"])
        :param standardize_st: phonopy structure, necessary for reading IR (Boolean)
        """
        for site_property in orig_st.site_properties:
            orig_st.remove_site_property(site_property)
        if vacuum_thickness:
            self.orig_st = Tools.set_vacuum(orig_st, vacuum_thickness)
        else:
            self.orig_st = orig_st
        self.defect_type = defect_type
        self.natom = natom
        self.vacuum_thickness = vacuum_thickness
        self.distort = None
        self.site_info = None
        self.defects = ChargedDefectsStructures(self.orig_st, cellmax=natom, antisites_flag=True).defects
        self.bulk_st = self.defects["bulk"]["supercell"]["structure"]
        self.NN_for_sudo_bulk = []

        if (defect_type[0] == "substitutions") or (defect_type[0] == "vacancies"):
            self.defect_entry = self.defects[self.defect_type[0]][self.defect_type[1]]
            self.defect_st = self.defect_entry["supercell"]["structure"].get_sorted_structure()
            self.defect_site_in_bulk = self.defect_entry["bulk_supercell_site"]
            self.defect_site_in_bulk_index = None
            self.NN = None

            if defect_type[0] == "substitutions":
                try:
                    self.defect_site_in_bulk_index = self.defect_st.index(self.defect_site_in_bulk)
                except ValueError:
                    self.defect_site_in_bulk_index = self.defect_st.index(self.defect_site_in_bulk.to_unit_cell())
                self.NN = [self.defect_st.index(self.defect_st[nn['site_index']])
                           for nn in CrystalNN().get_nn_info(self.defect_st, self.defect_site_in_bulk_index)]
                self.pmg_obj = Substitution(self.orig_st, self.defect_entry["unique_site"])

            elif defect_type[0] == "vacancies":
                try:
                    self.defect_site_in_bulk_index = self.bulk_st.index(self.defect_site_in_bulk)
                except ValueError:
                    self.defect_site_in_bulk_index = self.bulk_st.index(self.defect_site_in_bulk.to_unit_cell())
                self.NN = [self.defect_st.index(self.bulk_st[nn['site_index']])
                           for nn in CrystalNN().get_nn_info(self.bulk_st, self.defect_site_in_bulk_index)]
                self.pmg_obj = Vacancy(self.orig_st, self.defect_entry["unique_site"])

            self.nn_dist = dict(before=None, after=None)
            self.nn_dist["before"] = dict(zip([str(idx) for idx in self.NN], range(len(self.NN))))
            print("defect site coord. = orig: {} unit: {}".format(
                self.defect_site_in_bulk, self.defect_site_in_bulk.to_unit_cell()))

            self.defect_entry["supercell"].pop("structure")
            self.defect_entry["supercell"]["bulk"] = self.bulk_st


        elif defect_type[0] == "bulk":
            self.defect_entry = self.defects[self.defect_type[0]]
            self.defect_entry["supercell"].pop("structure")
            self.defect_st = None

        else:
            print("!!!Please insert substitutions, vacancies, or bulk!!!")

        if standardize_st and self.defect_st:
            self.defect_st, self.site_info = Tools.phonopy_structure(self.defect_st)

        if defect_type[0] == "substitutions":
            self.substitutions(distort, sub_on_side)

        elif defect_type[0] == "vacancies":
            self.vacancies(distort, sub_on_side)

    def substitutions(self, distort, substitution):
        bond_length = [self.defect_st.get_distance(self.defect_site_in_bulk_index, NN_index)
                       for NN_index in self.NN]
        bond_length = np.array(bond_length).round(3)

        self.nn_dist["before"] = dict(zip([str(idx) for idx in self.NN], bond_length))

        if substitution:
            self.make_complex(substitution)

        self.NN.append(self.defect_site_in_bulk_index)
        self.nn_dist["before"][str(self.defect_site_in_bulk_index)] = 0
        print("==" * 50, "\nBefore distortion: {}".format(self.nn_dist["before"]))

        if distort:
            self.move_sites(distort)
            bond_length = [self.defect_st.get_distance(self.defect_site_in_bulk_index, NN_index)
                           for NN_index in self.NN]
            bond_length = np.array(bond_length).round(3)
            self.nn_dist["after"] = dict(zip([str(idx) for idx in self.NN], bond_length))
            print("After distortion: {}\n{}".format(self.nn_dist["after"], "==" * 50))

    def vacancies(self, distort, substitution):
        bond_length = [self.bulk_st.get_distance(self.defect_site_in_bulk_index, NN_index['site_index'])
                       for NN_index in CrystalNN().get_nn_info(self.bulk_st, self.defect_site_in_bulk_index)]
        bond_length = np.array(bond_length).round(3)

        self.nn_dist["before"] = dict(zip([str(idx) for idx in self.NN], bond_length))

        if substitution:
            self.make_complex(substitution)

        print("==" * 50, "\nBefore distortion: {}".format(self.nn_dist["before"]))

        if distort:
            sudo_bulk = self.move_sites(distort)
            bond_length = [sudo_bulk.get_distance(self.defect_site_in_bulk_index, NN_index['site_index'])
                           for NN_index in CrystalNN().get_nn_info(sudo_bulk, self.defect_site_in_bulk_index)]
            bond_length = np.array(bond_length).round(3)
            self.nn_dist["after"] = dict(zip([str(idx) for idx in self.NN], bond_length))
            print("After distortion: {}\n{}".format(self.nn_dist["after"], "==" * 50))

    def move_sites(self, distort):
        self.distort = distort
        sudo_bulk = self.bulk_st.copy()
        if self.NN_for_sudo_bulk:
            NN_tot = zip(self.NN, self.NN_for_sudo_bulk)
        else:
            NN_tot = zip(self.NN, self.NN)
        for site, sudo_bulk_site in NN_tot:
            perturb = Tools.get_rand_vec(distort)
            self.defect_st.translate_sites([site], perturb, frac_coords=False)
            sudo_bulk.translate_sites([sudo_bulk_site], perturb, frac_coords=False)
        return sudo_bulk

    def move_origin_to_defect(self):
        # center_site_idx = self.NN[nn_idx]
        # center_site_coords = self.defect_st[center_site_idx].coords
        center_site_coords = self.defect_site_in_bulk.coords
        self.defect_st.translate_sites(range(self.defect_st.num_sites), -1*center_site_coords, frac_coords=False)

    def make_complex(self, substitution):
        self.defect_entry["complex"] = {"site": [], "site_specie": []}
        for sub, idx in zip(substitution, range(len(substitution))):
            self.defect_st.replace(self.NN[idx], sub)
            print("substitution coord = {}".format(self.NN[idx]))

            self.NN_for_sudo_bulk.append(self.NN[idx])

            self.defect_entry["complex"]["site"].append(self.defect_st[self.NN[idx]])
            self.defect_entry["complex"]["site_specie"].append(sub)

        self.defect_entry["defect_type"] = "complex"

        defect_sites_in_bulk = [self.defect_st[nn] for nn in self.NN]

        self.defect_st.sort()
        if self.defect_type[0] == "substitutions":
            self.defect_site_in_bulk_index = self.defect_st.index(self.defect_site_in_bulk)

        self.NN = [self.defect_st.index(nn) for nn in defect_sites_in_bulk]
        self.nn_dist["before"] = dict(zip([str(idx) for idx in self.NN], self.nn_dist["before"].values()))

if __name__ == '__main__':
    pc = Structure.from_file("example/H-BN.vasp")
    scaling = Tools.find_scaling_for_2d_defect(pc, 15)[0]
    area = scaling[0]*scaling[1]
    na = area*pc.num_sites
    thick = [20]
    distort = 0

    antisite_st = GenDefect(
        orig_st=pc,
        defect_type=("vacancies", 0),
        natom=na,
        vacuum_thickness=thick,
        sub_on_side=None,
        distort=0,
        standardize_st=False
    )

    vacancy = antisite_st.defect_st
    defect_entry_from_pycdt = antisite_st.defect_entry