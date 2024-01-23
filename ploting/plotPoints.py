import pandas as pd
import matplotlib.pyplot as plt
import sys


def plot_points(filename):
    # Load data from the CSV file
    data = pd.read_csv(filename, index_col=0)

    # Extract time values from the header
    with open(filename, "r") as file:
        header = file.readline().strip()
        # Remove the trailing comma and split
        pts = [float(p) for p in header[1:].split(",") if p]

    # Plot the data for each time with decreasing transparency
    for i, p in enumerate(pts):
        # Calculate alpha based on the index
        alpha = 1.0 - i / len(pts)
        plt.plot(data.index, data.iloc[:, i], label=f"Lx={p}", linewidth=2, alpha=alpha)
    plt.xlabel("t (s)")
    plt.ylabel("Temp (°C)")
    plt.title(r"Solutions at Different Points")
    plt.legend()
    plt.grid(False)
    plt.show()


if __name__ == "__main__":
    # Check if a filename is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python plotTimes.py <filename>")
        sys.exit(1)

    # Get the filename from the command-line argument
    filename = sys.argv[1]

    # Call the function to plot the data
    plot_points(filename)
