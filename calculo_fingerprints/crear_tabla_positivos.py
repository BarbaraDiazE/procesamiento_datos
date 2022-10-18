import pandas as pd
import os
familias_0 = [
    "LFA / ICAM",
    "Bromodominios / Histonas",
    "BCL2-Like / BAX",
    "XIAP / Smac",
    "LEDGF / IN",
    "CD80 / CD28",
    "MDM2-Like / P53",
    "CD4 / gp120"
]
familias= [
    "LFA / ICAM",
    "Bromodominios / Histonas",
    "BCL2-Like / BAX",
    "XIAP / Smac",
    "LEDGF / IN",
    "CD80 / CD28",
    "MDM2-Like / P53",
    "CD4 / gp120"
]
if __name__ == "__main__":

    df = dict()
    for i in familias_0:
        print(i)
        familias.remove(i)
        str1 = ""
        df[i] = [str1.join(familias)]
        familias = [
            "LFA / ICAM",
            "Bromodominios / Histonas",
            "BCL2-Like / BAX",
            "XIAP / Smac",
            "LEDGF / IN",
            "CD80 / CD28",
            "MDM2-Like / P53",
            "CD4 / gp120"
        ]
    print(df['CD80 / CD28'])
    for key, value in df.items():
        print("## ", key, " ###", value)
    df = pd.DataFrame.from_dict(df,
                                # columns = ["positivo", "negativo"],
                                orient="index"
                                )
    print(os.getcwd())
    print(df.iloc[0])
    df.to_csv("/home/babs/Desktop/positivos.tsv", sep="\t")
