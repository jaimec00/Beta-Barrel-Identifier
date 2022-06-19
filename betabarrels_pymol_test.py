
def color_beta():
    with open("pymol_barrel_test.txt", "r") as p:
        lines = p.readlines()
        pdbs = lines[0].split(",")
        color_commands = lines[1:]

        # Load the PDBs
        for pdb in pdbs:
            cmd.fetch(pdb)

        # Color everything green
        cmd.color("green")

        # Color sheets blue and barrels red
        for command in color_commands:
            command = command.split(', ')
            cmd.color(command[0], command[1])

cmd.extend("color_beta", color_beta)

color_beta()
