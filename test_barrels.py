# import pymol

def color_beta():

    with open("sheets.txt", "r") as s:
        lines = s.readlines()
        if len(lines) > 0:
            pdbs = lines[0].split(",")
            color_commands = lines[1:]
            for pdb in pdbs:
                cmd.fetch(pdb)
            # Color everything green
            cmd.color("green")
            # Color sheets blue
            for command in color_commands:
                cmd.color("blue", command)

    # Color barrels red
    with open("barrels.txt", "r") as b:
        color_commands = b.readlines()
        if len(color_commands) > 0:
            for command in color_commands:
                cmd.color("red", command)

cmd.extend("color_beta", color_beta)
