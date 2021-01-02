#! usr/bin/env python
import sys
import re

def getLength(full):
    count = 0
    if type(full) != list:
        Regex = r"\s"
        Replace = r""
        full = re.sub(Regex, Replace, full)
    for letter in full:
        count += 1
    return count

def checkArgs():
    if (getLength(sys.argv) != 6 or (getLength(sys.argv) == 6 and sys.argv[1] == "-l")) and (getLength(sys.argv) != 8 or sys.argv[6] != "-c"):
        print("Error: USAGE: python myprogram.py: option -i input.fasta -o output.txt")
        print("\tAvailable options:\n\t\t-g GC percent\n\t\t-r reverse complement")
        print("\t\t-s transcription\n\t\t-l translation\n\t\t-n count nucleotides")
        return
    if sys.argv[1] != "-g" and sys.argv[1] != '-s' and sys.argv[1] != '-l' and sys.argv[1] != '-r' and sys.argv[1] != '-n':
        print("Error: USAGE: python myprogram.py: option -i input.fasta -o output.txt")
        print("\tAvailable options:\n\t\t-g GC percent\n\t\t-r reverse complement")
        print("\t\t-s transcription\n\t\t-l translation\n\t\t-n count nucleotides")
        return
    if (getLength(sys.argv) == 8 and sys.argv[1] != "-l"):
        print("Error: USAGE: python myprogram.py: option -i input.fasta -o output.txt")
        print("\tAvailable options:\n\t\t-g GC percent\n\t\t-r reverse complement")
        print("\t\t-s transcription\n\t\t-l translation\n\t\t-n count nucleotides")
        return
    if (getLength(sys.argv) == 6 and (sys.argv[2] != "-i" or sys.argv[4] != "-o")):
        print("Error: USAGE: python myprogram.py: option -i input.fasta -o output.txt")
        print("\tAvailable options:\n\t\t-g GC percent\n\t\t-r reverse complement")
        print("\t\t-s transcription\n\t\t-l translation\n\t\t-n count nucleotides")
        return
    with open(sys.argv[3]) as inFile, open(sys.argv[5], "w") as outFile:
        count = 0
        if (sys.argv[1] == "-g"):
            header = ""
            seq = ""
            outFile.write("ID\tGC%\n")
            for line in inFile:
                line.strip()
                if count % 2 == 0:
                    header = remove_newLine(line) + "\t"
                    header = header[1:]
                    count += 1
                    continue
                else:
                    seq = prepare_seq(line)
                    count += 1
                outFile.write(calcGC(header, seq))
                seq = ""
            return
        elif (sys.argv[1] == "-r"):
            header = ""
            seq = ""
            for line in inFile:
                line.strip()
                if count % 2 == 0:
                    header = line
                    count += 1
                    continue
                else:
                    seq = prepare_seq(line)
                    count += 1
                outFile.write(revComp(header, seq))
                seq = ""
            return
        elif (sys.argv[1] == "-s"):
            header = ""
            seq = ""
            for line in inFile:
                line.strip()
                if count % 2 == 0:
                    header = line
                    count += 1
                    continue
                else:
                    seq = prepare_seq(line)
                    count += 1
                outFile.write(transcribe(header, seq))
                seq = ""
            return
        elif (sys.argv[1] == "-l"):
            header = ""
            seq = ""
            for line in inFile:
                line.strip()
                if count % 2 == 0:
                    header = line
                    count += 1
                    continue
                else:
                    seq = prepare_seq(line)
                    count += 1
                outFile.write(translate(header, seq, getCodons()))
                seq = ""
            return
        elif (sys.argv[1] == "-n"):
            header = ""
            seq = ""
            outFile.write("ID\tLength\tA(%A)\tC(%C)\tG(%G)\tT(%T)\n")
            for line in inFile:
                line.strip()
                if count % 2 == 0:
                    header = remove_newLine(line) + "\t"
                    header = header[1:]
                    count += 1
                    continue
                else:
                    seq = prepare_seq(line)
                    count += 1
                outFile.write(countNucs(header, seq))
                seq = ""
            return
        else:
            print("Error: USAGE: python myprogram.py: option -i input.fasta -o output.txt")
            print("\tAvailable options:\n\t\t-g GC percent\n\t\t-r reverse complement")
            print("\t\t-s transcription\n\t\t-l translation\n\t\t-n count nucleotides")
            return

def getCodons():
    retDict = {}
    with open(sys.argv[7]) as codonFile:
        for line in codonFile:
            line.strip()
            Regex = r"\n"
            Replace = r""
            line = re.sub(Regex, Replace, line)
            fields = line.split("\t")
            retDict[fields[0]] = fields[1]
        return retDict

def getCount(full, part):
    count = 0
    for letter in full:
        if part == letter:
            count += 1
    return count

def prepare_seq(line):
    retString = make_upper_case(line)
    Regex = r"\s"
    Replace = r""
    retString = re.sub(Regex, Replace, retString)
    return retString

def remove_newLine(line):
    Regex = r"\n"
    Replace = r""
    retstring = re.sub(Regex, Replace, line)
    return retstring

def make_upper_case(fullString):
    retString = ""
    for letter in fullString:
        if letter == 'a':
            retString += "A"
            continue
        elif letter == 'b':
            retString += "B"
            continue
        elif letter == 'd':
            retString += "D"
            continue
        elif letter  == 'e':
            retString += "E"
            continue
        elif letter == 'f':
            retString += "F"
            continue
        elif letter ==  'h':
            retString += "H"
            continue
        elif letter == 'i':
            retString += "I"
            continue
        elif letter == 'j':
            retString += "J"
            continue
        elif letter == 'k':
            retString += "K"
            continue
        elif letter == 'l':
            retString += "L"
            continue
        elif letter == 'm':
            retString += "M"
            continue
        elif letter == 'n':
            retString += "N"
            continue
        elif letter == 'o':
            retString += "O"
            continue
        elif letter == 'p':
            retString += "P"
        elif letter == 'q':
            retString += "Q"
            continue
        elif letter == 'r':
            retString += "R"
            continue
        elif letter == 's':
            retString += "S"
            continue
        elif letter == 'v':
            retString += "V"
            continue
        elif letter == 'w':
            retString += "W"
            continue
        elif letter == 'x':
            retString += "X"
            continue
        elif letter == 'y':
            retString += "Y"
            continue
        elif letter == 'z':
            retString += "Z"
            continue
        elif letter == 't':
            retString += "T"
            continue
        elif letter == 'g':
            retString += "G"
            continue
        elif letter == 'c':
            retString += "C"
            continue
        elif letter == 'u':
            retString += "U"
            continue
        else:
            retString += letter
            continue
    return retString



def convert_to_rna(value):
    text = make_upper_case(value)
    Regex = r"T"
    Replace = r"U"
    text = re.sub(Regex, Replace, text)
    Regex = r" "
    Replace = r""
    text = re.sub(Regex, Replace, text)
    return text

def validate_nucleotide_count(fullString):
    Regex = r" "
    Replace = r""
    fullString = re.sub(Regex, Replace, fullString)
    fullString = make_upper_case(fullString)
    if 'B' in fullString:
        return False
    elif 'D' in fullString:
        return False
    elif 'E' in fullString:
        return False
    elif 'F' in fullString:
        return False
    elif 'H' in fullString:
        return False
    elif 'I' in fullString:
        return False
    elif 'J' in fullString:
        return False
    elif 'K' in fullString:
        return False
    elif 'L' in fullString:
        return False
    elif 'M' in fullString:
        return False
    elif 'N' in fullString:
        return False
    elif 'O' in fullString:
        return False
    elif 'P' in fullString:
        return False
    elif 'Q' in fullString:
        return False
    elif 'R' in fullString:
        return False
    elif 'S' in fullString:
        return False
    elif 'U' in fullString:
        return False
    elif 'V' in fullString:
        return False
    elif 'W' in fullString:
        return False
    elif 'X' in fullString:
        return False
    elif 'Y' in fullString:
        return False
    elif 'Z' in fullString:
        return False
    elif '.' in fullString:
        return False
    elif '!' in fullString:
        return False
    elif ',' in fullString:
        return False
    elif '?' in fullString:
        return False
    elif '"' in fullString:
        return False
    elif '\'' in fullString:
        return False
    elif ':' in fullString:
        return False
    elif ';' in fullString:
        return False
    elif '0' in fullString:
        return False
    elif '1' in fullString:
        return False
    elif '2' in fullString:
        return False
    elif '3' in fullString:
        return False
    elif '4' in fullString:
        return False
    elif '5' in fullString:
        return False
    elif '6' in fullString:
        return False
    elif '7' in fullString:
        return False
    elif '8' in fullString:
        return False
    elif '9' in fullString:
        return False
    elif '@' in fullString:
        return False
    elif '#' in fullString:
        return False
    elif '$' in fullString:
        return False
    elif '%' in fullString:
        return False
    elif '^' in fullString:
        return False
    elif '&' in fullString:
        return False
    elif '*' in fullString:
        return False
    elif '(' in fullString:
        return False
    elif ')' in fullString:
        return False
    elif '-' in fullString:
        return False
    elif '_' in fullString:
        return False
    elif '+' in fullString:
        return False
    elif '=' in fullString:
        return False
    elif '{' in fullString:
        return False
    elif '}' in fullString:
        return False
    elif '[' in fullString:
        return False
    elif ']' in fullString:
        return False
    elif '|' in fullString:
        return False
    elif '"' in fullString:
        return False
    elif '\'' in fullString:
        return False
    elif '<' in fullString:
        return False
    elif '>' in fullString:
        return False
    elif '/' in fullString:
        return False
    return True

def validate_sequence(value):
    if not remove_newLine(value).isalpha():
        return False
    if value[0:3] != "ATG":
        return False
    if "TGA" not in value or (value.find("TGA") == value.rindex("TGA") and value.rindex("TGA") % 3 != 0):
        if "TAG" not in value or (value.find("TAG") == value.rindex("TAG") and value.rindex("TAG") % 3 != 0):
            if "TAA" not in value or (value.find("TAA") == value.rindex("TAA") and value.rindex("TAA") % 3 != 0):
                return False
    if 'B' in value:
        return False
    if 'D' in value:
        return False
    if 'E' in value:
        return False
    if 'F' in value:
        return False
    if 'H' in value:
        return False
    if 'I' in value:
        return False
    if 'J' in value:
        return False
    if 'K' in value:
        return False
    if 'L' in value:
        return False
    if 'M' in value:
        return False
    if 'N' in value:
        return False
    if 'O' in value:
        return False
    if 'P' in value:
        return False
    if 'Q' in value:
        return False
    if 'R' in value:
        return False
    if 'S' in value:
        return False
    if 'U' in value:
        return False
    if 'V' in value:
        return False
    if 'W' in value:
        return False
    if 'X' in value:
        return False
    if 'Y' in value:
        return False
    if 'Z' in value:
        return False
    if '.' in value:
        return False
    if '!' in value:
        return False
    if ',' in value:
        return False
    if '?' in value:
        return False
    if '"' in value:
        return False
    if '\'' in value:
        return False
    if ':' in value:
        return False
    if ';' in value:
        return False
    if '0' in value:
        return False
    if '1' in value:
        return False
    if '2' in value:
        return False
    if '3' in value:
        return False
    if '4' in value:
        return False
    if '5' in value:
        return False
    if '6' in value:
        return False
    if '7' in value:
        return False
    if '8' in value:
        return False
    if '9' in value:
        return False
    if '@' in value:
        return False
    if '#' in value:
        return False
    if '$' in value:
        return False
    if '%' in value:
        return False
    if '^' in value:
        return False
    if '&' in value:
        return False
    if '*' in value:
        return False
    if '(' in value:
        return False
    if ')' in value:
        return False
    if '-' in value:
        return False
    if '_' in value:
        return False
    if '+' in value:
        return False
    if '=' in value:
        return False
    if '{' in value:
        return False
    if '}' in value:
        return False
    if '[' in value:
        return False
    if ']' in value:
        return False
    if '|' in value:
        return False
    if '\'' in value:
        return False
    if '<' in value:
        return False
    if '>' in value:
        return False
    if '/' in value:
        return False
    return True

def translate(header, valid_sequence, codons_dictionary):
    retString = header
    if validate_sequence(valid_sequence):
        rnaSeq = convert_to_rna(valid_sequence)
        while getLength(rnaSeq) > 3:
            if rnaSeq[0:3] in codons_dictionary:
                if codons_dictionary[rnaSeq[0:3]] == "*":
                    break
                else:
                    retString += codons_dictionary[rnaSeq[0:3]]
                    rnaSeq = rnaSeq[3:]
            else:
                rnaSeq = rnaSeq[3:]
        return retString + "\n"
    else:
        return retString + "ERROR\n"

def countNucs(header, valid_sequence):
    fullString = header
    if validate_nucleotide_count(valid_sequence):
        fullCount = getLength(valid_sequence)
        aCount = getCount(valid_sequence, "A")
        tCount = getCount(valid_sequence, "T")
        cCount = getCount(valid_sequence, "C")
        gCount = getCount(valid_sequence, "G")
        aPercent = str((aCount / fullCount) * 100)
        tPercent = str((tCount / fullCount) * 100)
        cPercent = str((cCount / fullCount) * 100)
        gPercent = str((gCount / fullCount) * 100)
        aCount = str(aCount)
        tCount = str(tCount)
        cCount = str(cCount)
        gCount = str(gCount)
        fullCount = str(fullCount)
        fullString += fullCount + "\t" + aCount + "(" + aPercent + "%)\t" + cCount + "(" + cPercent + "%)\t" + gCount + "(" + gPercent + "%)\t" + tCount + "(" + tPercent + "%)\n"
    else:
        fullString += "ERROR\n"
    return fullString

def revComp(header, valid_sequence):
    if (validate_nucleotide_count(valid_sequence)):
        #if (valid_sequence == "TTTTA"):
         #   return header + "TAAAA\n"
        strand = header
        text = valid_sequence
        Regex = r"A"
        Replace = r"U"
        text = re.sub(Regex, Replace, text)
        Regex = r"C"
        Replace = r"J"
        text = re.sub(Regex, Replace, text)
        Regex = r"T"
        Replace = r"A"
        text = re.sub(Regex, Replace, text)
        Regex = r"U"
        Replace = r"T"
        text = re.sub(Regex, Replace, text)
        Regex = r"G"
        Replace = r"C"
        text = re.sub(Regex, Replace, text)
        Regex = r"J"
        Replace = r"G"
        text = re.sub(Regex, Replace, text)
        text = text[::-1]
        strand = strand + text + "\n"
        return strand
    else:
        return header + "ERROR\n"

def calcGC(header, valid_sequence):
    if (validate_nucleotide_count(valid_sequence)):
        retString = header
        fullCount = getLength(valid_sequence)
        count1 = getCount(valid_sequence, "C") + getCount(valid_sequence, "G")
        gcContent = str((count1 / fullCount) * 100)
        retString += gcContent + "%\n"
        return retString
    else:
        return header + "ERROR\n"

def transcribe(header, valid_sequence):
    retValue = header
    if (validate_nucleotide_count(valid_sequence)):
        for letter in valid_sequence:
            if letter == "T":
                retValue += "U"
            else:
                retValue += letter
        return retValue + "\n"
    else:
        return retValue + "ERROR\n"

checkArgs()