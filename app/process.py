import os
import sys
import time
sys.path.append(os.path.join(os.getcwd(), 'data'))

from Bio import SeqIO
from app.pysnpprocess import MpvProcess

import pandas as pd
import app.models as dbs

def raw_data_handler(raw, ref, gff):
    print("raw_handle")
    cwd = os.getcwd()
    updir = os.path.join(cwd, 'app', 'upload')

    uniqued = os.path.join(updir, 'uniqued.fasta')
    SeqIO.write(raw, uniqued, 'fasta')

    snps = os.path.join(updir, 'nucmer.snps')
    os.chdir(updir)
    os.system("nucmer --forward -p nucmer "+ ref + " "+ uniqued)
    os.system("show-snps nucmer.delta -T -l > "+ snps)

    process = MpvProcess(gff, ref)
    process.run(snps)

    os.remove(os.path.join(updir, 'uniqued.fasta'))
    os.chdir(cwd)

def insertion_clearer():
    updir = os.path.join(os.getcwd(), 'app', 'upload')
    os.remove(os.path.join(updir, 'nucmer.snps'))
    os.remove(os.path.join(updir, 'nucmer.delta'))

def check_lineage(lineage, species):
    local_lineage = lineage.unique().tolist()
    for i in local_lineage:
        obj, created = dbs.Lineage.objects.get_or_create(lineage = i)
        if(created):
            obj.species = species
            obj.save()

def sample_insertion(meta, species):
    seg = { 'seg' : False }
    print("sample_insert")

    if(meta.endswith('.xls') | meta.endswith('.xlsx')):
        mt = pd.read_excel(meta)
    else:
        mt = pd.read_csv(meta)

    if 'segment' in mt.columns:
        seg["seg"] = True
        seg["hash"] = {}
        for index, row in mt.iterrows():
            try:
                seg["hash"][row['segment']].append(row['ID'])
            except KeyError:
                seg["hash"].update({ row['segment']:[] })
                seg["hash"][row['segment']].append(row['ID'])

    meta_list = []
    check_lineage(mt['lineage'], species)
    mt['country'] = mt['country'].str.replace('Taiwan', 'China(Taiwan)')
    mt['country'] = mt['country'].str.replace('Hong Kong', 'China(Hong Kong)')
    mt['country'] = mt['country'].str.replace('Macau', 'China(Macau)')
    mt['country'] = mt['country'].str.replace('Viet Nam', 'Vietnam')
    mt['country'] = mt['country'].str.replace('USSR', 'Russia')
    mt['country'] = mt['country'].str.replace('Korea', 'South Korea')
    mt['country'] = mt['country'].str.replace('South South Korea', 'South Korea')
    mt['lineage'] = mt['lineage'].str.replace('Mixed', 'mixed')
    mt['lineage'] = mt['lineage'].str.replace('n1', 'N1')

    cdict = {}
    ldict = {}
    for i in list(dbs.Country.objects.all()):
        cdict[i.country] = i

    for i in list(dbs.Lineage.objects.all()):
        ldict[i.country] = i

    if(seg["seg"]):
        for i in mt.index:
            meta_list.append(dbs.Sample(id = mt['ID'][i], segment = mt['segment'][i], sample = mt['sample'][i], time = mt['time'][i], country = cdict[mt['country'][i]], lineage = ldict[mt['lineage'][i]], species = species))
    else:
        for i in mt.index:
            meta_list.append(dbs.Sample(id = mt['ID'][i], sample = mt['sample'][i], time = mt['time'][i], country = cdict[mt['country'][i]], lineage = ldict[mt['lineage'][i]], species = species))

    dbs.Sample.objects.bulk_create(meta_list, ignore_conflicts = True)
    return seg

def data_insertion(species, dataframe):
    print("data_insert")

    df = dataframe

    df['M_type'] = df['refvar'] + '->' + df['qvar']
    df['PM_type'] = df['refpos'].apply(str) + ':' + df['M_type']
    df['pro_variant'] = df['protein'] + ':' + df['variant']
    df = df.fillna(value = {'varclass' : 'NoArg', 'protein' : 'NA'})

    pos_list = []
    protein_index = {}
    sample_index = {}
    ex_list = []
    try:
        pkv = dbs.Position.objects.latest('pk').pk
    except dbs.Position.DoesNotExist:
        pkv = 0

    for i in species.protein_set.all():
        protein_index[i.protein] = i

    for i in dbs.Sample.objects.filter(id__in = df['sample'].tolist()):
        sample_index[i.id] = i

    count = 0
    for j in df.index:
        pkv = pkv + 1
        pos_list.append(dbs.Position(id = pkv, rpos = df['refpos'][j], rvar = df['refvar'][j], qpos = df['qpos'][j], qvar = df['qvar'][j], sample = sample_index[df['sample'][j]], protein = protein_index[df['protein'][j]]))
        ex_list.append(dbs.Extra(id = pkv, M_type = df['M_type'][j], PM_type = df['PM_type'][j], pro_variant = df['pro_variant'][j], variant = df['variant'][j], varclass = dbs.Extra.Varclass[df['varclass'][j]], position_id = pkv))
        count = count + 1
        if(count % 100000 == 0):
            dbs.Position.objects.bulk_create(pos_list)
            dbs.Extra.objects.bulk_create(ex_list)
            pos_list = []
            ex_list = []
            print(count / total)
    dbs.Position.objects.bulk_create(pos_list)
    dbs.Extra.objects.bulk_create(ex_list)
    print("inserted")

def data_entry(species):
    s = dbs.Species.objects.get(species = species)
    updir = os.path.join(os.getcwd(), 'app', 'upload')
    datadir = os.path.join(os.getcwd(), 'app', 'data')
    seg = False
    meta = ""
    raw = ""
    raws = []
    seen = set()

    for fd in os.listdir(updir):
        if(len(meta) != 0 & len(raw) != 0):
            break

        if(len(meta) == 0):
            if(fd.endswith('.csv')):
                meta = os.path.join(updir, fd)
                continue
            elif (fd.endswith('.xlsx')):
                meta = os.path.join(updir, fd)
                continue

        if(len(raw) == 0 & fd.endswith('.fasta')):
            raw = os.path.join(updir, fd)
            continue
    
    for record in SeqIO.parse(raw, 'fasta'):
        if str(record.id) not in seen:
            seen.add(str(record.id))
            raws.append(record)

    seg = sample_insertion(meta = meta, species = s)
    opt = ''
    outputfd = os.path.join(os.getcwd(), 'app', 'upload', 'output.csv')
    if(seg["seg"]):
        prefix = s.refseq + "_"
        output = []
        for key in seg["hash"]:
            templist = []
            for it in raws:
                if(it.id in seg["hash"][key]):
                    templist.append(it)
            
            raw_data_handler(raw = templist, ref = os.path.join(datadir, s.refseq, prefix + str(key) + '.fasta'), gff = os.path.join(datadir, s.refseq, prefix + str(key) + '.tsv'))
            output.append(pd.read_csv(outputfd))
            os.remove(outputfd)
        opt = pd.concat(output, ignore_index = True)
    else:
        raw_data_handler(raw = raws, ref = os.path.join(datadir, s.refseq), gff = os.path.join(datadir, s.annotation))
        opt = pd.read_csv(outputfd)
        os.remove(outputfd)
    
    data_insertion(species = s, dataframe = opt)
    insertion_clearer()
    #os.remove(meta)
    #os.remove(raw)
