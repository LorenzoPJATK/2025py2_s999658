#!/usr/bin/env python3
import sys,pandas as pd,matplotlib.pyplot as plt
from Bio import Entrez,SeqIO

class R:
 def __init__(s,e,k):Entrez.email,Entrez.api_key,Entrez.tool=e,k,'GenBankRetrieverFixed'
 def s(s,t):
  print(f"Searching for records for taxID: {t}")
  try:
   h=Entrez.efetch(db="taxonomy",id=t,retmode="xml")
   o=Entrez.read(h)[0]["ScientificName"]
   print(f"Organism found: {o} (TaxID: {t})")
   q=f"txid{t}[Organism]"
   h=Entrez.esearch(db="nucleotide",term=q,usehistory="y")
   r=Entrez.read(h)
   c=int(r["Count"])
   if not c:return print("No records found.")
   print(f"Records found: {c}")
   s.webenv,s.query_key,s.count=r["WebEnv"],r["QueryKey"],c
   return c
  except Exception as e:print(f"Error during search: {e}")
 def f(s,start=0,maxr=500):
  if not all(hasattr(s,a) for a in['webenv','query_key']):
   print("Run search_taxid() first.");return []
  try:
   h=Entrez.efetch(db="nucleotide",rettype="gb",retmode="text",retstart=start,retmax=maxr,webenv=s.webenv,query_key=s.query_key)
   return list(SeqIO.parse(h,"gb"))
  except Exception as e:print(f"Error during fetch: {e}");return []

def gen_csv(rs,f):
 df=pd.DataFrame([{"Accession":r.id,"Length":len(r.seq),"Description":r.description} for r in rs])
 df.to_csv(f,index=False);print(f"CSV report saved to: {f}");return df

def gen_plot(df,f):
 df=df.sort_values(by="Length",ascending=False).reset_index(drop=True)
 plt.figure(figsize=(12,6))
 plt.plot(df["Accession"],df["Length"],marker='o')
 plt.xticks(rotation=90,fontsize=8)
 plt.xlabel("GenBank Accession")
 plt.ylabel("Sequence Length")
 plt.title("Sequences sorted by length")
 plt.tight_layout();plt.savefig(f)
 print(f"Chart saved as: {f}")

def main():
 e="s999658@pjwstk.edu.pl";k="1350188c32df2b0d825bd26bc9ada74a5a09"
 t="9606"
 minl = int(input("Enter minimum sequence length: "))
 maxl = int(input("Enter maximum sequence length: "))
 maxr = 100
 r=R(e,k);c=r.s(t)
 if not c:sys.exit("No records found. Program terminated.")
 print(f"\nDownloading up to {maxr} records...")
 rs=r.f(0,maxr)
 fr=[x for x in rs if minl<=len(x.seq)<=maxl]
 print(f"Records filtered by length [{minl}, {maxl}]: {len(fr)}")
 if not fr:sys.exit("No records meet the length criteria. Program terminated.")
 csv,f=f"taxid_{t}_report.csv",f"taxid_{t}_chart.png"
 df=gen_csv(fr,csv);gen_plot(df,f)

if __name__=="__main__":main()