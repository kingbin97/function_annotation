import sys

file_query = sys.argv[1]
file_database = sys.argv[2]
file_out = sys.argv[3]

h = {}

with open(file_database, 'r') as database:
    for line in database:
        s, d = line.strip().split('\t')[0], line.strip().split('\t')[1]
        h[s] = d

with open(file_out, 'w') as out:
    with open(file_query, 'r') as query:
        for line in query:
            line_parts = line.strip().split('\t')
            q = line_parts[1]
            if q in h:
                GO_ID = h[q].replace(" ", "")  # 去掉空格
                out.write(f"{line_parts[0]}\t{q}\t{GO_ID}\n")
            else:
                pass #out.write(line)  # 如果没有匹配到，不输出|直接写入原始行内容
