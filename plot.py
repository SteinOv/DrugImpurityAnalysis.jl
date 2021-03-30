from sys import argv, exit
import csv

def main():

    # Data filename optional
    if len(argv) == 1:
        data = "234MMC der 1.csv"
    elif len(argv) == 2:
        data = argv[1]
    else:
        print("Usage: python3 plot.py {data}")
        exit(1)







    # Print intensities of chosen masses
    with open(f"data/{data}", 'r') as file:
        data = csv.reader(file)
        masses = [82, 182, 303]
        mass_col = [i - 41 + 4 for i in masses]

        next(data)

        intensity_list = []
        for i in range(len(masses)):
            intensity_list.append(0)

        for line in data:
            for i, mass in enumerate(mass_col):
                intensity_list[i] += int(float(line[mass]))

        print(intensity_list)

if __name__ == "__main__":
    main()