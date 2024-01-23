import pandas as pd
import matplotlib.pyplot as plt
import sys


def plot_times(filename):
    # Load data from the CSV file
    data = pd.read_csv(filename, index_col=0)

    # Extract time values from the header
    with open(filename, "r") as file:
        header = file.readline().strip()
        # Remove the trailing comma and split
        times = [float(t) for t in header[1:].split(",") if t]

    # Plot the data for each time with decreasing transparency
    for i, time in enumerate(times):
        # Calculate alpha based on the index
        alpha = 1.0 - i / len(times)
        plt.plot(
            data.index, data.iloc[:, i], label=f"t={time}s", linewidth=2, alpha=alpha
        )

    # Customize the plot
    plt.xlabel("x (m)")
    plt.ylabel("Temp (°C)")
    plt.title("Solutions at Different Times")
    plt.grid(False)
    plt.legend(ncol=3, fontsize=8)
    plt.show()


if __name__ == "__main__":
    # Check if a filename is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python plotTimes.py <filename>")
        sys.exit(1)

    # Get the filename from the command-line argument
    filename = sys.argv[1]

    # Call the function to plot the data
    plot_times(filename)
