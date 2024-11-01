import pickle
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix
from imblearn.over_sampling import SMOTE



# create a custom colourmap for use throughout
colors = ["#5e8fb4", "#FFFFFF", "#e6b4b0"]
n_bins = 500  # Discretizes the interpolation into bins
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=n_bins)


def load_dataframe(inputfile):

    # with open(snakemake.input.embedding_df, "rb") as input_file:
    #     embedding_df = pickle.load(input_file)

    # rb is read bytes
    # this loads it into a pandas df
    with open(inputfile, "rb") as input_file:
        embedding_df = pickle.load(input_file)

    embedding_df['Clade'] = embedding_df['info'].apply(tag_node)

    return embedding_df


def tag_node(info):
    with open('./data/reportdata/combined_nodes/NR1_ids.txt', 'r') as file:
        nr1_names = set(line.strip() for line in file)

    with open('./data/reportdata/combined_nodes/NR4_ids.txt', 'r') as file:
        nr4_names = set(line.strip() for line in file)

    if info in nr1_names:
        return 'NR1'
    elif info in nr4_names:
        return 'NR4'
    else:
        return 'Other'
    


embeddings_df = load_dataframe('./data/ancestor_embedding_combined_df.csv')
embeddings_df['Clade'] = embeddings_df['info'].apply(tag_node)
df = embeddings_df.drop(columns=['info', 'sequence', 'model_name', 
                                 'protbert_max_embedding', 'protbert_mean_embedding', 
                                 'protbert_weighted_embedding'])

df = df.join(pd.DataFrame(df.pop('protbert_cls_embedding').tolist(), index=df.index))


df.to_csv('logregtest_df.csv')


print(df)
print()
print()


print((df['Clade'] == 'NR1').sum())
print((df['Clade'] == 'NR4').sum())


# Step 1: Encode the Clade column if necessary
le = LabelEncoder()
df['Clade_encoded'] = le.fit_transform(df['Clade'])

# Step 2: Define features (X) and target (y)
X = df.drop(columns=['Clade', 'Clade_encoded'])  # Drop non-numerical columns
y = df['Clade_encoded']

# Step 3: Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# apply smote to correct class imbalance
smote = SMOTE(random_state=42)
X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)


# Step 4: Scale the features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train_resampled)
X_test_scaled = scaler.transform(X_test)

# Step 6: Train the logistic regression model
model = LogisticRegression(max_iter=1000)  # You can increase max_iter if necessary
model.fit(X_train_scaled, y_train_resampled)

# Step 7: Predict and evaluate the model
y_pred = model.predict(X_test_scaled)

# Evaluate accuracy and confusion matrix
accuracy = accuracy_score(y_test, y_pred)
conf_matrix = confusion_matrix(y_test, y_pred)

print(f"Accuracy: {accuracy}")
print(f"Confusion Matrix: \n{conf_matrix}")