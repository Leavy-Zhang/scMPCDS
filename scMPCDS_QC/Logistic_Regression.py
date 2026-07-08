import os,re
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score
from sklearn.feature_selection import SelectKBest, chi2, RFECV
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns


np.random.seed(42)


def simulate_snp_allele_fraction(n_samples, n_snps, n_causal_snps=5):

    
    snp_data = np.random.beta(2, 3, size=(n_samples, n_snps))  
    
    
    for i in range(n_causal_snps):
        
        case_mask = np.random.binomial(1, 0.5, size=n_samples).astype(bool)
        snp_data[case_mask, i] = np.random.beta(3, 2, size=np.sum(case_mask))  
    
    
    
    causal_effects = np.random.normal(2, 0.5, size=n_causal_snps)
    
    
    risk_score = np.dot(snp_data[:, :n_causal_snps], causal_effects)
    
    
    prob = 1 / (1 + np.exp(-(risk_score - 5)))
    
    
    phenotype = np.random.binomial(1, prob)
    
    
    while len(np.unique(phenotype)) < 2:
        phenotype = np.random.binomial(1, prob)
    
    
    snp_columns = [f'SNP_{i+1}' for i in range(n_snps)]
    df = pd.DataFrame(snp_data, columns=snp_columns)
    df['phenotype'] = phenotype
    
    return df


def load_data(file_path=None, nonzero_frac = 0.01):
    """
    load data
    file_path: if None, generate simulated training data
    """
    if file_path:
        df = pd.read_csv(file_path, header=0 ,index_col=0)
        frac = (df.iloc[:,:(df.shape[1]-1)] != 0).sum()/(df.shape[1]-1)
        frac = (df.iloc[:,:(df.shape[1]-1)] != 0).sum()/(df.shape[1]-1)
        df = pd.concat([df[frac.index[frac >= nonzero_frac]], df.iloc[:, -1]], axis =1)
    else:
        
        n_samples = 1000
        n_snps = 100
        n_causal_snps = 10
        df = simulate_snp_allele_fraction(n_samples, n_snps, n_causal_snps)
        print(f"generated {n_samples} samples, ratios of of {n_snps} SNP sites")
        df.to_csv('snp_allele_fraction_matrix.csv', index=False)
        print("mutation frequency matrix of SNP is saved as snp_allele_fraction_matrix.csv")
    
    return df


def preprocess_data(df):
    """
    data preprocessing
    """
    
    X = df.iloc[:,:-1]
    y = df.iloc[:, -1]
    
    return X, y


def select_features(X, y, method='chi2', k=20):
    """
    selected key features
    method: 'chi2' or 'rfe'
    k: number of features returned
    """
    if method == 'chi2':
        selector = SelectKBest(chi2, k=k)
        X_selected = selector.fit_transform(X, y)
        selected_features = X.columns[selector.get_support()]
    elif method == 'rfe':
        estimator = LogisticRegression(max_iter=1000)
        selector = RFECV(estimator, step=1, cv=5, scoring='accuracy')
        X_selected = selector.fit_transform(X, y)
        selected_features = X.columns[selector.get_support()]
    else:
        raise ValueError("Method must be 'chi2' or 'rfe'")
    
    return X_selected, selected_features


def build_logistic_regression(X, y, prefix = "NA"):
    """
    build and train model
    """
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    
    param_grid = {
        'C': [0.001, 0.01, 0.1, 1, 10, 100],
        'penalty': ['l1', 'l2'],
        'solver': ['liblinear']
    }
    
    grid_search = GridSearchCV(LogisticRegression(max_iter=1000), param_grid, cv=5, scoring='accuracy')
    grid_search.fit(X_train_scaled, y_train)
    
    
    best_model = grid_search.best_estimator_
    print(f"optimal super parameter: {grid_search.best_params_}")
    
    
    y_pred = best_model.predict(X_test_scaled)
    y_pred_proba = best_model.predict_proba(X_test_scaled)[:, 1]
    
    
    accuracy = accuracy_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_pred_proba)
    
    print("Model estimation result:")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"AUC: {auc:.4f}")
    print("Classification report:")
    print(classification_report(y_test, y_pred))
    a = classification_report(y_test, y_pred, output_dict=True)
    a = pd.DataFrame(a).transpose()
    a.to_csv("path/to/" + prefix+".fusion_matrix.csv")
    return  best_model, accuracy, auc


def identify_key_snps(model, features):
    """
    Key SNP screening using model coefficients
    """
    coefficients = model.coef_[0]
    snp_importance = pd.DataFrame({
        'SNP': features,
        'Coefficient': coefficients,
        'Absolute_Coefficient': np.abs(coefficients)
    })
    
    
    snp_importance = snp_importance.sort_values('Absolute_Coefficient', ascending=False)
    
    return snp_importance


def optimize_model(X, y, key_snps, top_n=10):
    """
    Model Optimization using screened key SNPs
    top_n: number of top key SNP
    """
    
    top_snps = key_snps['SNP'].head(top_n).tolist()
    X_top = X[top_snps]
    
    
    print(f"using top {top_n} SNPs to rebuld models...")
    optimized_model, accuracy, auc = build_logistic_regression(X_top, y)
    
    return optimized_model, accuracy, auc, top_snps


def main():


    prefix = "mito_LR_analysis"
    df = load_data(file_path="path/to/"+prefix+".csv", nonzero_frac = 0.01)


    df.to_csv("path/to/"+ prefix+ '.filtered.csv', index=False)
    print("Mutation fraction matrix is saved as snp_allele_fraction_matrix.csv")


    print("Data preprocessing...")
    X, y = preprocess_data(df)


    print("key feature selection...")
    X_selected, selected_features = select_features(X, y, method='chi2', k=20)
    print(f"{len(selected_features)} key SNP features are selected")


    print("Build logistic regression model...")
    model, accuracy, auc = build_logistic_regression(X_selected, y)


    print("Screeing key SNP...")
    key_snps = identify_key_snps(model, selected_features)
    print("Information of top 10 key SNP:")
    print(key_snps.head(10))


    key_snps.to_csv("path/to/"+ prefix + '.key_snps_from_allele_fraction.csv', index=False)
    print("Screened key SNP is saved as key_snps_from_allele_fraction.csv")


    print("Optimizing model...")
    optimized_model, opt_accuracy, opt_auc, top_snps = optimize_model(X, y, key_snps, top_n=10)


    plt.figure(figsize=(12, 8))
    sns.barplot(x='SNP', y='Absolute_Coefficient', data=key_snps.head(20))
    plt.xticks(rotation=90)
    plt.title('Top 20 Key SNPs by Coefficient')
    plt.tight_layout()
    plt.savefig("path/to/"+ prefix + '.key_snps_allele_fraction_plot.pdf')
    print("Visualization of key SNP is save as key_snps_allele_fraction_plot.pdf")
    
    fusion_mat = pd.read_csv("path/to/" + prefix+".fusion_matrix.csv", header=0,index_col=0)
    tmp = fusion_mat.iloc[:2,:3]
    tmp.index = ['low','high']
    fig, ax= plt.subplots(figsize=(3,2))
    sns.heatmap(tmp, cmap=("Blues"),annot=True,ax=ax)
    ax.set_title(prefix)
    fig.savefig(os.path.join("path/to/", prefix+"fusion.pdf"))

if __name__ == "__main__":
    main()



