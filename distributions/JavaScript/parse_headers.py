import os

include_directory = os.fsencode("../.././include")
function_types = ["int", "double", "long", "void"]
function_names = []
inc_file_names = []
function_counter = 0


def browse_files(directory):
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        inc_file_names.append(filename)
        function_names.append(parse_file("../.././include/" + filename))


def parse_file(filename):
    result = []

    with open(filename, "r") as file:
        lines = file.readlines()

    for line in lines:
        split_line = line.split()
        if len(split_line) >= 2 and split_line[0] in function_types:
            if split_line[1][0] == "A":
                result.append(split_line[1][:-1])

    return result


def create_script(script_name):
    i = 0

    with open(script_name, "w") as file:
       
        file.write("\nfuns='['\n")

        for header in function_names:
            file.write("\n#----" + inc_file_names[i] + "\n")

            for fun in header:
                file.write("funs+='\"_" + fun + "\",' \n")

            i += 1

        file.write("\nfuns+=']'\n")


def print_missing_files(directory):
    print("Files not parsed: ")

    for file in os.listdir(directory):
        filename = os.fsdecode(file)

        if filename not in inc_file_names:
            print(filename)


def main():
    browse_files(include_directory)
    create_script("funs.txt")


if __name__ == "__main__":
    main()
