# load libraries
import pandas as pd 
import numpy as np 
import seaborn as sns
from fpdf import FPDF
import sys
import io
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

###########################################User Mods ######################################################################
# Prep df for analysis
# step 1: user variables
df = pd.read_excel('/Users/Shrini/Desktop/BCA_Assay/BCA_29_1_2_2025.xlsx', index_col = 0)
csv_output = '/Users/Shrini/Desktop/BCA_Assay/BCA_29_1_2_2025_result.csv'
pdf_output_path = '/Users/Shrini/Desktop/BCA_Assay/BCA_29_1_2_2025_report.pdf'
starting_dilution = 2
sample_no = 4 # var1: sample #

#### Verify these settings and change if needed
std_conc = [1, 0.5, 0.25, 0.125, 0.0625, 0.0313]
dilution_factor = [2, 4, 6, 8, 10, 12]
###################################################################################################################

sample_no = sample_no + 2 # accounting for the two standards
df = df.iloc[:, :sample_no] # subsetting to useful columns only

#normalize for blank
blank = round(df.loc["G"].mean(), 3) # generate blank and 

df = df.drop(index=["G", "H"], errors="ignore")# Drop rows G and H
df = df - blank # normalize by removing blank


df["std_conc"] = std_conc # add my standard concentration
##### Generate a function to plot trendlines

def calculate_trendline(X, Y):
    X_reshaped = X.values.reshape(-1, 1)  # Reshape for sklearn
    model = LinearRegression()
    model.fit(X_reshaped, Y)
    Y_pred = model.predict(X_reshaped)
    r2 = r2_score(Y, Y_pred)
    slope = model.coef_[0]
    intercept = model.intercept_
    return Y_pred, slope, intercept, r2

# Calculate trend line for col1
Y_pred1, slope1, intercept1, r2_1 = calculate_trendline(df['std_conc'], df[1])

# Calculate trend line for col2
Y_pred2, slope2, intercept2, r2_2 = calculate_trendline(df['std_conc'], df[2])

# Plotting
plt.figure(figsize=(8,6))
plt.scatter(df['std_conc'], df[1], label='standard_1 Data', color='blue', alpha=0.7)
plt.plot(df['std_conc'], Y_pred1, color='blue', linestyle='dashed',
         label=f'Std1 Trendline: y = {slope1:.4f}x + {intercept1:.4f}\n$R^2$ = {r2_1:.4f}')

plt.scatter(df['std_conc'], df[2], label='Standard2 Data', color='red', alpha=0.7)
plt.plot(df['std_conc'], Y_pred2, color='red', linestyle='dashed',
         label=f'std2 Trendline: y = {slope2:.4f}x + {intercept2:.4f}\n$R^2$ = {r2_2:.4f}')

# Labels and Title
plt.xlabel("std_conc (X-axis)")
plt.ylabel("Y values")
plt.title("Scatter Plot with Two Trendlines")
plt.legend()
plt.grid(True)

# Show plot
plt.show()

# compare the two slopes and select one formulae
if r2_1 > r2_2:
    best_slope, best_intercept = slope1, intercept1
    best_r2 = r2_1
    best_label = "col1"
else:
    best_slope, best_intercept = slope2, intercept2
    best_r2 = r2_2
    best_label = "col2"

print(f"The best fitting line is for {best_label}:")
print(f"y = {best_slope:.4f}x + {best_intercept:.4f}")
print(f"R^2 = {best_r2:.4f}")

# now back calculate and give the me concentrations for each 
y_columns = df.drop(columns=[1, 2, 'std_conc'])
solved_x = (y_columns - best_intercept) / best_slope
solved_x['dilution_factor'] = dilution_factor # now multiply with dilution factor
solved_x_scaled = solved_x.mul(solved_x['dilution_factor'], axis=0)
# Ensure the 'dilution_factor' column is not included in the result
solved_x_scaled = solved_x_scaled.drop(columns=['dilution_factor'], errors='ignore')
print("")
print("Concentration calculate for each sample in all wells")
print(solved_x_scaled)
print()
# select 3 values to compute the mean protein abundance of the three closest values

# Function to find the three most similar values in a column and compute their mean
def mean_of_three_most_similar(series):
    sorted_values = np.sort(series.values)  # Sort values in ascending order
    diffs = np.diff(sorted_values)  # Compute absolute differences between consecutive values
    min_diff_indices = np.argsort(diffs)[:2]  # Get indices of the two smallest differences

    # Select three values, ensuring we always get three
    first_index = min_diff_indices[0]  # Start from the first closest pair
    if first_index + 2 < len(sorted_values):  
        selected_values = sorted_values[first_index:first_index + 3]
    else:  
        selected_values = sorted_values[-3:]  # Pick the last three if at the end

    return np.mean(selected_values), selected_values

# Apply function to each column in solved_x_scaled
results = solved_x_scaled.apply(mean_of_three_most_similar, axis=0)

# Extract means and selected values separately
closest_values_mean = results.apply(lambda x: x[0])  # Extract mean values
selected_values = results.apply(lambda x: x[1])  # Extract selected values
selected_values

print("Based on our analyses these were the three closest results")
print(selected_values)
print()
print("")
print("")
print("Shown here are the concentration of protein in mg/ml")
print(closest_values_mean)
point_result = closest_values_mean


#################################################Slope analysis##################################################################
###################################################################################################################
###################################################################################################################

filtered_df = df.where(df >= blank, np.nan)
filtered_df = filtered_df.drop(columns=[1, 2, 'std_conc'], errors="ignore")


# Reassign std_conc while ensuring row alignment is maintained
# Use original index from df to avoid row loss

# Function to calculate slope using linear regression while ignoring NaNs
def calculate_slope_ignore_nans(X, Y):
    valid_mask = X.notna() & Y.notna()  # Mask to ignore NaNs in both X and Y
    X_valid = X[valid_mask]
    Y_valid = Y[valid_mask]

    if len(X_valid) < 3 or len(Y_valid) < 3:  # Ensure enough points for regression
        return np.nan

    X_reshaped = X_valid.values.reshape(-1, 1)  
    model = LinearRegression()
    model.fit(X_reshaped, Y_valid)
    return model.coef_[0]

filtered_df['std_conc'] = std_conc 

# Ensure std_conc column exists
if 'std_conc' not in filtered_df.columns:
    raise ValueError("Column 'std_conc' is required but missing in filtered_df")

# Calculate slopes while ensuring std_conc is properly filtered
slope_results = filtered_df.drop(columns=['std_conc'], errors="ignore").apply(
    lambda col: calculate_slope_ignore_nans(filtered_df['std_conc'], col), axis=0
)

# Normalize by best_slope
final_results = slope_results / best_slope

slope_results = final_results * starting_dilution
slope_results


slope_results["source"] = "slope_results"
point_result["source"] = "point_result"

# Merge on index
merged_result = pd.concat([slope_results, point_result], axis=1, join="outer")
print("This is the final result accounting for both")
print(merged_result)
merged_result.to_csv(csv_output)


###################################################################################################################
#############################################Report Generation #################################################################
###################################################################################################################

# Generate plot
plot_path = "/Users/Shrini/Desktop/SCTLD_Proteomics/tmp_plot.png"
plt.figure(figsize=(8,6))
plt.scatter(df['std_conc'], df[1], label='standard_1 Data', color='blue', alpha=0.7)
plt.plot(df['std_conc'], Y_pred1, color='blue', linestyle='dashed',
         label=f'Std1 Trendline: y = {slope1:.4f}x + {intercept1:.4f}\n$R^2$ = {r2_1:.4f}')

plt.scatter(df['std_conc'], df[2], label='Standard2 Data', color='red', alpha=0.7)
plt.plot(df['std_conc'], Y_pred2, color='red', linestyle='dashed',
         label=f'Std2 Trendline: y = {slope2:.4f}x + {intercept2:.4f}\n$R^2$ = {r2_2:.4f}')

plt.xlabel("std_conc (X-axis)")
plt.ylabel("Y values")
plt.title("Scatter Plot with Two Trendlines")
plt.legend()
plt.grid(True)
plt.savefig(plot_path)
plt.close()

# Capture print output
output_buffer = io.StringIO()
sys.stdout = output_buffer

# Rerun print statements
print(f"The best fitting line is for {best_label}:")
print(f"y = {best_slope:.4f}x + {best_intercept:.4f}")
print(f"R^2 = {best_r2:.4f}")
print("")
print("Concentration calculated for each sample in all wells")
print(solved_x_scaled)
print("")
print("Based on our analyses these were the three closest results")
print(selected_values)
print("")
print("Shown here are the concentration of protein in mg/ml for slope and point analysis")
print(merged_result)

# Restore standard output
sys.stdout = sys.__stdout__

# Get captured text
output_text = output_buffer.getvalue()

# Create PDF
pdf = FPDF()
pdf.set_auto_page_break(auto=True, margin=15)
pdf.add_page()
pdf.set_font("Arial", size=12)

# Add text to PDF
pdf.multi_cell(0, 10, output_text)

# Add the plot image
pdf.add_page()
pdf.image(plot_path, x=10, y=20, w=180)

# Save PDF
pdf.output(pdf_output_path)

# Return the PDF path
print(f"PDF saved at: {pdf_output_path}")

