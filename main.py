from enum import Enum
import pymol_scripting


class Protein(Enum):
    Phe = "Phe"
    Leu = "Leu"
    Ser = "Ser"
    Tyr = "Tyr"
    Stop = "Stop"
    Cys = "Cys"
    Trp = "Trp"
    Pro = "Pro"
    His = "His"
    Gln = "Gln"
    Arg = "Arg"
    Ile = "Ile"
    Start = "Start"
    Thr = "Thr"
    Asn = "Asn"
    Lys = "Lys"
    Val = "Val"
    Ala = "Ala"
    Asp = "Asp"
    Glu = "Glu"
    Gly = "Gly"


def Protein_to_Blast(p: Protein) -> str:
    d = {
        Protein.Phe: "F",
        Protein.Leu: "L",
        Protein.Ser: "S",
        Protein.Tyr: "Y",
        Protein.Stop: "*",
        Protein.Cys: "C",
        Protein.Trp: "W",
        Protein.Pro: "P",
        Protein.His: "H",
        Protein.Gln: "Q",
        Protein.Arg: "R",
        Protein.Ile: "I",
        Protein.Start: "M",
        Protein.Thr: "T",
        Protein.Asn: "N",
        Protein.Lys: "K",
        Protein.Val: "V",
        Protein.Ala: "A",
        Protein.Asp: "D",
        Protein.Glu: "E",
        Protein.Gly: "G"
    }
    return d[p]


def Resolve_Triple(s: str) -> Protein:
    if len(s) != 3: raise Exception()
    CodeSonne = {
        "uuu": Protein.Phe,
        "uuc": Protein.Phe,
        "uua": Protein.Leu,
        "uug": Protein.Leu,
        "ucu": Protein.Ser,
        "ucc": Protein.Ser,
        "uca": Protein.Ser,
        "ucg": Protein.Ser,
        "uau": Protein.Tyr,
        "uac": Protein.Tyr,
        "uaa": Protein.Stop,
        "uag": Protein.Stop,
        "ugu": Protein.Cys,
        "ugc": Protein.Cys,
        "uga": Protein.Stop,
        "ugg": Protein.Trp,
        "cuu": Protein.Leu,
        "cuc": Protein.Leu,
        "cua": Protein.Leu,
        "cug": Protein.Leu,
        "ccu": Protein.Pro,
        "ccc": Protein.Pro,
        "cca": Protein.Pro,
        "ccg": Protein.Pro,
        "cau": Protein.His,
        "cac": Protein.His,
        "caa": Protein.Gln,
        "cag": Protein.Gln,
        "cgu": Protein.Arg,
        "cgc": Protein.Arg,
        "cga": Protein.Arg,
        "cgg": Protein.Arg,
        "auu": Protein.Ile,
        "auc": Protein.Ile,
        "aua": Protein.Ile,
        "aug": Protein.Start,
        "acu": Protein.Thr,
        "acc": Protein.Thr,
        "aca": Protein.Thr,
        "acg": Protein.Thr,
        "aau": Protein.Asn,
        "aac": Protein.Asn,
        "aaa": Protein.Lys,
        "aag": Protein.Lys,
        "agu": Protein.Ser,
        "agc": Protein.Ser,
        "aga": Protein.Arg,
        "agg": Protein.Arg,
        "guu": Protein.Val,
        "guc": Protein.Val,
        "gua": Protein.Val,
        "gug": Protein.Val,
        "gcu": Protein.Ala,
        "gcc": Protein.Ala,
        "gca": Protein.Ala,
        "gcg": Protein.Ala,
        "gau": Protein.Asp,
        "gac": Protein.Asp,
        "gaa": Protein.Glu,
        "gag": Protein.Glu,
        "ggu": Protein.Gly,
        "ggc": Protein.Gly,
        "gga": Protein.Gly,
        "ggg": Protein.Gly
    }
    return CodeSonne[s]


def load_file(filename):
    with open(filename, 'r') as f:
        data = f.readlines()[1:]
        data = ''.join(data).replace('\n', '').lower().replace("t", "u")
    return data


def Find_ORFS(filename):
    data = load_file(filename)
    orfs = []
    for i in range(3, len(data)+1):
        if Resolve_Triple(data[i - 3:i]) == Protein.Start:
            current = [Protein.Start]
            end = i - 3
            for j in range(i + 3, len(data)+1, 3):
                if Resolve_Triple(data[j - 3:j]) == Protein.Stop:
                    end = j - 3
                    break
                current.append(Resolve_Triple(data[j - 3:j]))

            else:
                #print("No stop codon found")
                continue
            if end - i < 20:
                continue
            orfs.append((current, i - 3, end))

    return [("".join(map(Protein_to_Blast, x[0])), x[1:]) for x in orfs]


def compare(orf1, orf2):
    orfres = []
    for i in range(len(orf1)):
        for j in range(len(orf2)):
            if orf1[i][0] in orf2[j][0]:
                if orf2[j] not in orfres:
                    orfres.append(orf2[j])
                break
            elif orf2[j][0] in orf1[i][0]:
                if orf1[i] not in orfres:
                    orfres.append(orf1[i])
                break

    print(len(orf1), len(orf2), len(orfres))
    return orfres


def get_mutations(file1, file2):
    file1 = load_file(file1)
    file2 = load_file(file2)
    mutations = [] # position

    for i in range(len(file1)):
        if file1[i] != file2[i]:
            mutations.append((i, file1[i], file2[i]))

    return mutations

if __name__ == '__main__':
    print(get_mutations("mutations/ORF15_original.fa", "mutations/ORF15_mutations.fa"))
    orf1 = Find_ORFS("mutations/ORF15_original.fa")
    orf2 = Find_ORFS("mutations/ORF15_mutations.fa")
    print(orf1)
    print(orf2)


"""
if __name__ == '__main__':
    print(Find_ORFS("input/1980.fasta"))
    exit()
    res = compare(Find_ORFS("input/1980.fasta"), Find_ORFS("input/2024.fasta"))

    print("\n".join(map(str, res)))
"""
gut = ["MRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQGEPTAHCCPWPPPATPCSWRSHPAWAEGGRRLPPSRGSGALF",
       #  "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQGEPTAHCCPWPPPATPCSWRSHPAWAEGGRRLPPSRGSGALF"
       "MGSETIKPAGAQQPSALQDRLHQKRPSSRSVPRAFASGGLRVPGWLDPRPQLCSREDVAGLVKHVGVSPGAPRQGTWPSACLSPACLPDHCPSAMALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQGEPTAHCCPWPPPATPCSWRSHPAWAEGGRRLPPSRGSGALF"]

# 'MGSETIKPAGAQQPSALQDRLHQKRPSSRSVPRAFASGGLRVPGWLDPRPQLCSREDVAGLVKHVGVSPGAPRQGTWPSACLSPACLPDHCPSAMALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQGEPTAHCCPWPPPATPCSWRSHPAWAEGGRRLPPSRGSGALF', (2141, 2735))
# ('MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQGEPTAHCCPWPPPATPCSWRSHPAWAEGGRRLPPSRGSGALF', (238, 550))
##('MGRRGQEAATQQGVRCTFLKRSSLGHVLKVTSSLWPSQNLSLRTVLASAAPRYIRGWARSSLHSPLKQMPRSPFLHPHLMTADSSVLLSKVLGDLGSQGAPRCLPLGEHPITPGGGRGCLPEWARPLSPASRQLHSQEMGKMLGTGPGEKYWDHLFRLPL', (2680, 3160))
# ('MPRSPFLHPHLMTADSSVLLSKVLGDLGSQGAPRCLPLGEHPITPGGGRGCLPEWARPLSPASRQLHSQEMGKMLGTGPGEKYWDHLFRLPL', (2884, 3160))
# ('MTADSSVLLSKVLGDLGSQGAPRCLPLGEHPITPGGGRGCLPEWARPLSPASRQLHSQEMGKMLGTGPGEKYWDHLFRLPL', (2917, 3160))
# ('MGKMLGTGPGEKYWDHLFRLPL', (910, 976)) - vlt gut
# ('MWALGPVGPHPVWVTLPLTWVQPGWRWVGVRPRAGGQAGTVSP', (1011, 1140))
