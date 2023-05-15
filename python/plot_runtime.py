import matplotlib.pyplot as plt
import numpy as np

# Data
c_max = [10, 15, 20, 25, 30, 35, 40]
kop1 = [0.035,
        0.078,
        0.12,
        0.17,
        0.2,
        0.28,
        0.35]
kop3 = [0.41,
        0.64,
        1.02,
        1.41,
        1.75,
        1.87,
        3.39]
kop6 = [2.19,
        3.03,
        4.49,
        5.33,
        9.02,
        12.65,
        15.05]

# Plotting
plt.figure()

# Plotting KOP-1_OUR_P
plt.plot(c_max, kop1, marker='o', linestyle=':', label='KOP_B-1')

# Plotting KOP-3_OUR_P
plt.plot(c_max, kop3, marker='o', linestyle='--', label='KOP_B-3')

plt.plot(c_max, kop6, marker='o', label='KOP_B-6')

# Set logarithmic scale for the y-axis
plt.yscale('log')
plt.legend(loc='upper left')
# Set the x and y axis labels
plt.xlabel('C_max in seconds')
plt.ylabel('Avg. runtime in seconds (log)')
for tick in plt.gca().get_yticks():
    plt.axhline(tick, color='gray', linestyle=':', linewidth=0.5)

# Set the y-axis limits
plt.ylim(0.01, 10000)

# Add a legend
plt.legend()

# Show the plot
plt.show()
