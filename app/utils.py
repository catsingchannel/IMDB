from django.db.models import Count
from django.db import connection
from numpy import unique
from app.models import *
from collections import ChainMap, Counter
from Bio import SeqIO
from app.primer_design import *
import time, math, itertools
import pandas as pd
import os
import json

def sample_filter(req):
    q = req
    q["species"] = Species.objects.get(species = q["species"]).id

    if('all' != q["segment"]):
        if not q["segment"]:
            q["segment"] = 0
        else:
            q["segment"] = int(q["segment"])
        raw = Sample.objects.filter(species__id = q["species"], segment = q["segment"], time__range = q["date"])
    else:
        raw = Sample.objects.filter(species__id = q["species"], time__range = q["date"])

    if('all' not in q["lineage"]):
        raw = raw.filter(lineage__in = Lineage.objects.filter(lineage__in = q["lineage"]))

    if('all' not in q["locate"]):
        regionlist = ['Asia', 'Europe', 'North America', 'South America', 'Africa', 'Oceania']
        countrylist = []
        if(bool(set(regionlist) & set(q["locate"]))):
            for it in q["locate"]:
                if it in regionlist:
                    for itt in Country.objects.filter(region = Country.Regions[it.replace(" ", "").upper()]).values("country"):
                        countrylist.append(itt['country'])
                else:
                    countrylist.append(it)
            countrylist = list(set(countrylist))
        else:
            countrylist = q["locate"]
        raw = raw.filter(country__in = Country.objects.filter(country__in = countrylist))

    return raw

def country_ridge_plot(req):
    raw = sample_filter(req)

    first = False
    data = {}
    data["columns"] = list(raw.values('country').annotate(total = Count('country')).order_by('-total')[:10].values_list('country', flat = True))
    data["data"] = []

    for it in data["columns"]:
        df = list(Position.objects.filter(sample__in = raw.filter(country_id = it), rpos__range = req['genomerange']).values_list('qpos', flat = True))
        if (first) :
            for itt in range(len(df)):
                data['data'][itt][it] = df[itt]
        else:
            for itt in df:
                temp = {}
                temp[it] = itt
                data['data'].append(temp)
            first = True

    return data

def mutation_classes(req):
    raw = sample_filter(req)

    data = {}
    bucket = {}

    start = time.time()
    pset = Position.objects.filter(sample__in = raw)
    #ids = list(pset.values_list('id', flat=True))
    mpslist = list(pset.values_list('sample', flat = True))
    end = time.time()
    print(end - start)
    dvalues, mpscounts = unique(mpslist, return_counts = True)
    mpsvalues, mpscounts = unique(mpscounts.tolist(), return_counts = True)
    mpssort = sorted(dict(zip(mpsvalues.tolist(), mpscounts.tolist())).items())

    nelist = []
    tnelist = list(pset.values_list('rvar', 'rpos', 'qvar'))
    for it in tnelist:
        nelist.append(it[0] + str(it[1]) + it[2])
    nevalues, necounts = unique(nelist, return_counts = True)
    nesort = sorted(dict(zip(nevalues.tolist(), necounts.tolist())).items(), key=lambda x:x[1], reverse=True)

    data['NucleoEvents'] = []
    for idx in nesort[:10]:
        data['NucleoEvents'].append({ 'Mutation' : idx[0], 'Occ' : idx[1] })

    data['MutPerSample'] = []
    for idx in mpssort:
        data['MutPerSample'].append({ 'Count' : str(idx[0]), 'Occ' : idx[1] })

    start = time.time()
    dset = list(Extra.objects.filter(position__in = pset).values('varclass', 'M_type'))
    end = time.time()
    print(end - start)
    vclist, mtlist = zip(*map(dict.values, dset))
    pvlist = Extra.objects.filter(position__in = pset).exclude(varclass__in = ['SS', 'EX']).values_list('pro_variant', flat=True)
    vcvalues, vccounts = unique(vclist, return_counts = True)
    mtvalues, mtcounts = unique(mtlist, return_counts = True)
    pvvalues, pvcounts = unique(pvlist, return_counts = True)
    vcsort = sorted(dict(zip(vcvalues.tolist(), vccounts.tolist())).items(), key=lambda x:x[1], reverse=True)
    mtsort = sorted(dict(zip(mtvalues.tolist(), mtcounts.tolist())).items(), key=lambda x:x[1], reverse=True)
    pvsort = sorted(dict(zip(pvvalues.tolist(), pvcounts.tolist())).items(), key=lambda x:x[1], reverse=True)

    data['VarClasses'] = []
    vcdict = dict(Extra.Varclass.choices)
    for idx in vcsort[:10]:
        data['VarClasses'].append({ 'Varclass' : vcdict[idx[0]], 'Occ' : idx[1] })

    data['VarType'] = []
    for idx in mtsort[:10]:
        data['VarType'].append({ 'Type' : idx[0], 'Occ' : idx[1] })

    data['ProEvents'] = []
    for idx in pvsort[:10]:
        data['ProEvents'].append({ 'Mutation' : idx[0], 'Occ' : idx[1] })
    return data

def single_gene(req, raw):
    data = {}
    data['total'] = raw.count()
    mtype_num = int(req["mtypeNum"])

    if(req['custom']):
        pset = list(Position.objects.filter(sample__in = raw, rpos__range = req['genomerange']).values('rvar', 'rpos', 'qvar'))
    else:
        p = Protein.objects.get(species__id = req['species'], protein = req["protein"])
        pset = list(Position.objects.filter(sample__in = raw, protein = p).values('rvar', 'rpos', 'qvar'))

    rvar, rpos, qvar = zip(*map(dict.values, pset))
    glist = [i+":"+j+":"+k for i,j,k in zip(rvar, list(map(str, rpos)), qvar)]
    mlist = [i+"->"+j for i,j in zip(rvar, qvar)]
    gvalues, gcounts = unique(glist, return_counts = True)
    mvalues, mcounts = unique(mlist, return_counts = True)
    gsort = sorted(dict(zip(gvalues.tolist(), gcounts.tolist())).items(), key=lambda x:x[1], reverse=True)
    msort = sorted(dict(zip(mvalues.tolist(), mcounts.tolist())).items(), key=lambda x:x[1], reverse=True)[:mtype_num]


    data["proteinBar"] = []
    data["proteinCount"] = {}
    data["proteinCount"]["data"] = []
    msort = [i[0] for i in msort]
    data["proteinCount"]["typeList"] = msort
    count = 0
    for idx in gsort:
        substr = idx[0].split(':')
        if(count < 10):
            data["proteinBar"].append({ "Mutation" : substr[0] + substr[1] + substr[2], "Occ" : idx[1] })
            count += 1
        mtype = substr[0] + "->" + substr[2]
        if(mtype in msort):
            data["proteinCount"]["data"].append({ "M_type" : mtype, "Position" : int(substr[1]), "proOcc" : idx[1] })

    return data

def mutation_rate(req):
    raw = sample_filter(req)
    data = single_gene(req, raw)

    if(req["global"]):
        s = Species.objects.get(id = req["species"])
        p = Protein.objects.filter(species = s)
        pro = list(p.values('id', 'protein', 'length'))
        pset = Position.objects.filter(sample__in = raw, protein__in = p)
        lenhash = { d["protein"]:d["length"] for d in pro }
        phash = { d["id"]:d["protein"] for d in pro }
        pcount = { phash[d["protein"]]:d["total"] for d in list(pset.values('protein').annotate(total=Count('sample', distinct=True))) }

        try:
            discard = [ pcount.pop(key) for key in ["3'UTR", "5'UTR", "intergenic"] ]
        except Exception as e:
            pass

        data["mutationRate"] = []
        for key, value in pcount.items():
            data["mutationRate"].append({ "Protein" : key, "MutationRate" : (value / lenhash[key]) })

    return data

def global_profile(req):
    data = {}
    ptype_num = int(req["ptypeNum"])
    raw = sample_filter(req)
    data['total'] = raw.count()

    p = Protein.objects.get(species__id = req['species'], protein = req["protein"])
    pset = Position.objects.filter(sample__in = raw, protein = p)
    eset = Extra.objects.filter(position__in = pset).exclude(varclass__in = ['SS', 'EX']).values('pro_variant', 'variant', 'PM_type')
    pv, v, glist = zip(*map(dict.values, list(eset)))
    pvc = [ i.split(":")[0] + '/' + j for i, j in zip(glist, pv) ]

    pvalues, pcounts = unique(pv, return_counts = True)
    cvalues, ccounts = unique(pvc, return_counts = True)
    vvalues, vcounts = unique(v, return_counts = True)
    csort = sorted(dict(zip(cvalues.tolist(), ccounts.tolist())).items(), key=lambda x:x[1], reverse=True)
    psort = sorted(dict(zip(pvalues.tolist(), pcounts.tolist())).items(), key=lambda x:x[1], reverse=True)[:ptype_num]
    vsort = sorted(dict(zip(vvalues.tolist(), vcounts.tolist())).items(), key=lambda x:x[1], reverse=True)

    data['proEvents'] = []
    for idx in vsort[:10]:
        data['proEvents'].append({ 'Mutation' : idx[0], 'Occ' : idx[1] })

    data["proteinProfile"] = {}
    data["proteinProfile"]["data"] = []
    psort = [i[0] for i in psort]
    data["proteinProfile"]["typeList"] = psort
    for idx in csort:
        substr = idx[0].split('/')
        ptype = substr[1]
        if(ptype in psort):
            data["proteinProfile"]["data"].append({ "pro_variant" : ptype, "Position" : int(substr[0]), "proOcc" : idx[1] })

    return data

def upset_venn(req):
    data = {}
    data["venn"] = {}
    data["venn"]["data"] = []
    data["venn"]["labels"] = []
    data["upset"] = {}
    data["upset"]["combinations"] = []
    data["upset"]["sets"] = []
    ctype_num = int(req["ctypeNum"])
    raw = sample_filter(req)
    req['muts'] = set(req['muts'])

    start = time.time()
    p = Protein.objects.get(species__id = req['species'], protein = req["protein"])
    pset = list(Position.objects.filter(sample__in = raw, protein = p).values('id', 'sample'))

    pid, sid = zip(*map(dict.values, pset))
    end = time.time()
    print(end - start)
    start = time.time()
    eset = list(Extra.objects.filter(id__in = pid).exclude(varclass__in = ["EX","SS"]).values('id', 'variant'))
    eid, variant = zip(*map(dict.values, eset))
    end = time.time()
    print(end - start)

    start = time.time()
    vvalues, vcounts = unique(variant, return_counts = True)
    vsort = sorted(dict(zip(vvalues.tolist(), vcounts.tolist())).items(), key=lambda x:x[1], reverse=True)[:ctype_num]
    end = time.time()
    print(end - start)
    for id, size in vsort:
        data["upset"]["sets"].append({ 'setId' : id, 'size' : size })
    vsort = [i[0] for i in vsort]

    svdict = {}
    sdict = { k:v for k,v in zip(pid, sid) }
    for id, value in zip(eid, variant):
        if(value in vsort):
            idx = sdict[id]
            if(idx not in svdict):
                svdict[idx] = set()
            svdict[idx].add(value)

    if(len(req['muts']) > 1):
        pset = list(Position.objects.filter(sample__in = raw).values('id', 'sample'))
        pid, sid = zip(*map(dict.values, pset))
        eset = list(Extra.objects.filter(id__in = pid, pro_variant__in = req['muts']).values('id', 'pro_variant'))
        if(len(eset) > 0):
            eid, variant = zip(*map(dict.values, eset))
        else:
            eid = ()
            variant = ()
        eset = list(Extra.objects.filter(id__in = pid, PM_type__in = req['muts']).values('id', 'PM_type'))
        if(len(eset) > 0):
            peid, mutation = zip(*map(dict.values, eset))
        else:
            peid = ()
            mutation = ()
        eid = eid + peid
        variant = variant + mutation

        venn_svdict = {}
        sdict = { k:v for k,v in zip(pid, sid) }
        for id, value in zip(eid, variant):
            idx = sdict[id]
            if(idx not in venn_svdict):
                venn_svdict[idx] = set()
            venn_svdict[idx].add(value)

    cdict = {}
    for key, values in svdict.items():
        idx = frozenset(values)
        if(idx not in cdict):
            cdict[idx] = []
        cdict[idx].append(key)

    count = 0
    for key, value in cdict.items():
        data["upset"]["combinations"].append({ 'id' : count, 'sets' : list(key), 'samples' : value })
        count += 1

    if(len(req['muts']) > 1):
        cdict = {}
        sdict = {}
        r = (2 ** (len(req['muts']) + 1))
        for s in range(1, len(req['muts']) + 1):
            for combs in itertools.combinations(req['muts'], s):
                cdict[frozenset(combs)] = 0
                sdict[frozenset(combs)] = r
            r = r / 2

        for key, values in venn_svdict.items():
            idx = frozenset(values)
            if(idx not in cdict):
                pass
            else:
                cdict[idx] += 1 

        id = 0
        for key, value in cdict.items():
            if(len(key) == 1):
                label = str(next(iter(key))) + "\n" + str(value)
            else:
                label = str(value)

            data["venn"]["data"].append({"sets" : list(key), 'size' : sdict[key]})
            data["venn"]["labels"].append({"sets" : list(key), 'label' : label})


    data["upset"]["combinations"] = sorted(data["upset"]["combinations"], key = lambda x: len(x["samples"]), reverse = True)

    return data

def assays(req):
    tmp = ""
    profile = {}

    if(req["custom"]):
        tmp = req["template"]
    else:
        s = Species.objects.get(species = req["species"])
        datadir = os.path.join(os.getcwd(), 'app', 'data')
        refseq = ""
        annot = ""

        if(req["range"] == False):
            req["segment"] = Protein.objects.filter(protein = req["protein"], species = s)[0].segment

        if(req["segment"] != 0):
            prefix = s.refseq + "_"
            refseq = os.path.join(datadir, s.refseq, prefix + str(req["segment"]) + '.fasta')
            annot = os.path.join(datadir, s.refseq, prefix + str(req["segment"]) + '.tsv')
        else:
            refseq = os.path.join(datadir, s.refseq)
            annot = os.path.join(datadir, s.annotation)
        refseq = list(SeqIO.parse(refseq, 'fasta'))[0]
        refseq = str(refseq.seq)

        raw = sample_filter(req)
        pset = []
        if(req["range"]):
            tmp = refseq[(int(req["genomerange"][0]) - 1):int(req["genomerange"][1])]
            req['start'] = int(req["genomerange"][0]) - 1
            pset = list(Position.objects.filter(sample__in = raw, rpos__range = req['genomerange']).values("id", "sample__id", "rpos"))
        else:
            annot = pd.read_csv(annot, index_col = 0, sep='\t')
            annot.fillna("NA", inplace=True)
            annot = annot.loc[annot["gene"] == req["protein"]]
            tmp = refseq[(annot["start"].iloc[0] - 1):(annot["end"].iloc[0])]
            req['start'] = annot["start"].iloc[0] - 1
            req['start'] = int(req['start'])

            pset = list(Position.objects.filter(sample__in = raw, protein = Protein.objects.get(protein = req["protein"], species = s)).values("id", "sample__id", "rpos"))

        ps_dict = {}
        ep_dict = {}
        if(len(pset) > 0):
            pid, sid, rpos = zip(*map(dict.values, pset))
            ps_dict = { k:v for k, v in zip(pid, sid) }
            ep_dict = { k:v for k, v in zip(pid, rpos) }

        pset = []
        for key, value in ep_dict.items():
            try:
                profile[value].append(ps_dict[key])
            except KeyError:
                profile.update({ value:[] })
                profile[value].append(ps_dict[key])

        ps_dict = {}
        ep_dict = {}

    param = req
    param["seq"] = tmp
    param["profile"] = profile

    return primer_design(**param)
