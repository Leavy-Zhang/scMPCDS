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
        print(f"生成了 {n_samples} 个样本，{n_snps} 个SNP位点的突变比例数据")
        df.to_csv('snp_allele_fraction_matrix.csv', index=False)
        print("SNP突变频率矩阵已保存到 snp_allele_fraction_matrix.csv")
    
    return df


def preprocess_data(df):
    """
    数据预处理
    """
    
    X = df.iloc[:,:-1]
    y = df.iloc[:, -1]
    
    return X, y


def select_features(X, y, method='chi2', k=20):
    """
    selected key features
    method: 'chi2' 或 'rfe'
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
    print(f"最佳超参数: {grid_search.best_params_}")
    
    
    y_pred = best_model.predict(X_test_scaled)
    y_pred_proba = best_model.predict_proba(X_test_scaled)[:, 1]
    
    
    accuracy = accuracy_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_pred_proba)
    
    print("模型评估结果:")
    print(f"准确率: {accuracy:.4f}")
    print(f"AUC: {auc:.4f}")
    print("分类报告:")
    print(classification_report(y_test, y_pred))
    a = classification_report(y_test, y_pred, output_dict=True)
    a = pd.DataFrame(a).transpose()
    a.to_csv("path/to/" + prefix+".fusion_matrix.csv")
    return  best_model, accuracy, auc


def identify_key_snps(model, features):
    """
    基于模型系数筛选关键SNP
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
    使用关键SNP优化模型
    top_n: 选择前n个关键SNP
    """
    
    top_snps = key_snps['SNP'].head(top_n).tolist()
    X_top = X[top_snps]
    
    
    print(f"使用前 {top_n} 个关键SNP重新构建模型...")
    optimized_model, accuracy, auc = build_logistic_regression(X_top, y)
    
    return optimized_model, accuracy, auc, top_snps


def main():


    prefix = "mito_LR_analysis"
    df = load_data(file_path="path/to/"+prefix+".csv", nonzero_frac = 0.01)


    df.to_csv("path/to/"+ prefix+ '.filtered.csv', index=False)
    print("SNP突变频率矩阵已保存到 snp_allele_fraction_matrix.csv")


    print("数据预处理...")
    X, y = preprocess_data(df)


    print("进行特征选择...")
    X_selected, selected_features = select_features(X, y, method='chi2', k=20)
    print(f"选择了 {len(selected_features)} 个重要SNP特征")


    print("构建逻辑回归模型...")
    model, accuracy, auc = build_logistic_regression(X_selected, y)


    print("筛选关键SNP...")
    key_snps = identify_key_snps(model, selected_features)
    print("前10个关键SNP:")
    print(key_snps.head(10))


    key_snps.to_csv("path/to/"+ prefix + '.key_snps_from_allele_fraction.csv', index=False)
    print("关键SNP已保存到 key_snps_from_allele_fraction.csv")


    print("优化模型...")
    optimized_model, opt_accuracy, opt_auc, top_snps = optimize_model(X, y, key_snps, top_n=10)


    plt.figure(figsize=(12, 8))
    sns.barplot(x='SNP', y='Absolute_Coefficient', data=key_snps.head(20))
    plt.xticks(rotation=90)
    plt.title('Top 20 Key SNPs by Coefficient')
    plt.tight_layout()
    plt.savefig("path/to/"+ prefix + '.key_snps_allele_fraction_plot.pdf')
    print("关键SNP可视化已保存到 key_snps_allele_fraction_plot.pdf")
    
    fusion_mat = pd.read_csv("path/to/" + prefix+".fusion_matrix.csv", header=0,index_col=0)
    tmp = fusion_mat.iloc[:2,:3]
    tmp.index = ['low','high']
    fig, ax= plt.subplots(figsize=(3,2))
    sns.heatmap(tmp, cmap=("Blues"),annot=True,ax=ax)
    ax.set_title(prefix)
    fig.savefig(os.path.join("path/to/", prefix+"fusion.pdf"))

if __name__ == "__main__":
    main()



