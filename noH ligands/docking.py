from vina import Vina

v = Vina(sf_name='vina', cpu=0, seed=0)

v.set_receptor(rigid_pdbqt_filename='6lu7noC.pdbqt')
v.set_ligand_from_file('noH_000519.pdbqt')
v.compute_vina_maps(center=[-7, 12, 67.7], box_size=[22, 22, 22], spacing=1.0)

v.dock(exhaustiveness=32, n_poses=10)
v.write_poses('000519_docked.pdbqt',n_poses=10, overwrite=True)

