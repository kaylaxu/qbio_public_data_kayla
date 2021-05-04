import cptac  # import cptac to download cptac protein and clinical data
import pandas as pd  # import pandas for
import matplotlib.pyplot as plt  # import matplotlib for boxplot
import seaborn as sns  # import seaborn for boxplot
from scipy import stats
from statannot import add_stat_annotation

cptac.download(dataset="Brca")   # download breast cancer dataset
br = cptac.Brca()  # save data in br variable

protein_data = br.get_proteomics()  # save proteomic data
protein_data = protein_data.droplevel(1, axis=1)     # remove multi index
clinical_data = br.get_clinical()  # save clinical data

esr1 = protein_data["ESR1"]  # save ESR1 protein expression column
clinical_data["ER.IHC.Score"] = clinical_data["ER.IHC.Score"].fillna("Not reported")   # fill in null values

er_mask = clinical_data["ER.IHC.Score"] == "3+"  # 3+ is ER-positive, create mask of ER-positive patients
patients = esr1[er_mask]  # apply mask to protein expression data
ages = clinical_data["Age.in.Month"][er_mask]/12  # calculate ages in years
er_positive_patients = pd.DataFrame(patients, columns=["ESR1"])  # create new dataframe
er_positive_patients["Age"] = ages  # apply ages to new column in dataframe
category = []  # set categories list
for age in er_positive_patients["Age"]:  # for each age in the age column
    if age < 50:  # if it's less than 50, append young/patient is considered young
        category.append("Young")
    else:  # if more than/equal to 50, append old/patient is considered old
        category.append("Old")
category = pd.array(category)  # make category a pandas array
er_positive_patients["Age_category"] = category   # create new column in df titled age category
plt.figure()
df = er_positive_patients  # set df to data frame
x = "Age_category"   # set x to age_category
y = "ESR1"  # set y to esr1 counts
order = ["Young", "Old"]  # order of x axis
ax = sns.boxplot(data=df, x=x, y=y, order=order)  # create boxplot
add_stat_annotation(ax, data=df, x=x, y=y, order=order,
                                   box_pairs=[("Young", "Old")],
                                   test='Mann-Whitney', text_format='star',
                                   loc='outside', verbose=2)     # calculate p-value
plt.savefig("/PATHWAY/esr1_protein_expression_young_old_boxplot.png")

# calculate regression
slope, intercept, r_value, p_value, std_err = stats.linregress(er_positive_patients["Age"], er_positive_patients["ESR1"])
plt.figure()
# create scatterplot w/regression line
ax = sns.regplot(x="Age", y="ESR1", data=er_positive_patients, color='b', line_kws={'label': "y={0:.1f}x+{1:.1f}".format(slope, intercept)})
# plot legend
ax.legend()
plt.savefig("/PATHWAY/esr1_protein_expression_young_old_linear.png")
