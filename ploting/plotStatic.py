import pandas as pd
import matplotlib.pyplot as plt


def plot_data(filename):
    # Load data from the CSV file
    data = pd.read_csv(filename)

    # Extract x values from the DataFrame
    x_values = data["x"]

    # Plot the 'sol' columns
    plt.plot(x_values, data["sol"], "b", label="numeric", linewidth=2)
    plt.plot(x_values, data["solEx"], "r--", label="exact", linewidth=2)

    # Customize the plot
    plt.xlabel("x (m)")
    plt.ylabel("Temp (°C)")
    plt.title("Exact vs Numerical Solutions")
    plt.legend()
    plt.grid(False)

    # Show the plot
    plt.show()


if __name__ == "__main__":
    import sys

    # Check if a filename is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <filename>")
        sys.exit(1)

    # Get the filename from the command-line argument
    filename = sys.argv[1]

    # Call the function to plot the data
    plot_data(filename)
