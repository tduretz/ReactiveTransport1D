using DataFrames, CSV

# Read data into a dataframe
df_Eugster87 = (CSV.read("data/MatrixBruciteEugster87.csv", DataFrame))

@show df_Eugster87

# Convert to Matrix
Matrix(df_Eugster87[:,2:end])

# Extract labels
labels = names(df_Eugster87)