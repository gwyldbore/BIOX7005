import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder
import seq_utils
import pickle
import joblib
from imblearn.over_sampling import SMOTE



def load_dataframe(inputfile):
    with open(inputfile, "rb") as input_file:
        embedding_df = pickle.load(input_file)

    embedding_df['Clade'] = embedding_df['info'].apply(seq_utils.tag_node)

    return embedding_df

def main():    
    # Load the df with the ancestor sequences
    embedding_df = load_dataframe(snakemake.input.ancestor_embeddings)

    # drop the columns i don't want
    df = embedding_df.drop(columns=['info', 'sequence', 'model_name', 
                                'protbert_max_embedding', 'protbert_cls_embedding', 
                                'protbert_weighted_embedding'])
    
    # transform into each feature as a single column
    df = df.join(pd.DataFrame(df.pop('protbert_mean_embedding').tolist(), index=df.index))

    df.columns = df.columns.astype(str)
    print('dataframe columns:', df)


    # split features into x and y
    X_train = df.drop(columns=["Clade"])
    y_train = df["Clade"]  

    # encode the clade labels
    label_encoder = LabelEncoder()
    y_train_encoded = label_encoder.fit_transform(y_train)

    # Apply SMOTE to the training data to correct class imbalance
    smote = SMOTE(random_state=42)
    X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train_encoded)

    # Train logistic regression on the resampled data
    logistic_model = LogisticRegression(max_iter=1000)
    logistic_model.fit(X_train_resampled, y_train_resampled)

    # save the trained model for later use
    joblib.dump((logistic_model, label_encoder), snakemake.output.model_output)



if __name__ == "__main__":
    main()