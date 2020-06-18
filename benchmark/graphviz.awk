BEGIN {
    printf("digraph {\n");
    if (red_file != "") {
        while (getline<red_file) {
            gsub("[-/\\.|]", "_", $1);
            red_nodes[$1] = 1;
            #printf("%s[style=filled,color=red];\n", $1);
        }
    }
    if (red_node != "") {
            gsub("[-/\\.|]", "_", red_node);
            red_nodes[red_node] = 1;
            #printf("%s[style=filled,color=red];\n", red_node);
    }
    if (mapping_file != "") {
        while(getline<mapping_file) {
            gsub("[-/\\.|]", "_", $1);
            gsub("[-/\\.|]", "_", $2);
            if ($2 != -1) {
                #mapping_nodes[$1] = $2"\t"$3;
                mapping_nodes[$1"\t"$2] = $3;
            }
        }
    }
    if (ref_file != "") {
        i=0
        while(getline<ref_file) {
            if (i % 2 == 1) {
                ref_tbl[(i-1)/2] = $0;
            }
            ++i;
        }
    }
} {
    if ($1 == "VT") {
        gsub("[-/\\.|]", "_", $2);
        LEN_TBL[$2] = length($3);
        for (i in ref_tbl) {
            pos = index(ref_tbl[i], $3);
            while (pos != 0) {
                if (POS_TBL[$2] == "") {
                    POS_TBL[$2] = i"_"pos;
                } else {
                    POS_TBL[$2] = POS_TBL[$2]"_"i"_"pos;
                }
                x = index(substr(ref_tbl[i], pos + 1), $3);
                if (x == 0) {
                    break;
                }
                pos = pos + x;
            }
        }
        if ($2 in red_nodes) {
            printf("%s_%d_%s[style=filled,color=red];\n", $2, length($3), POS_TBL[$2]);
        } else {
            printf("%s_%d_%s;\n", $2, length($3), POS_TBL[$2]);
        }
    } else if ($1=="ED" && $5 - $4 + 1 >= minOverlap) {
        gsub("[-/\\.|]", "_", $2);
        gsub("[-/\\.|]", "_", $3);
        EDGE_NODE[$2] = 1;
        EDGE_NODE[$3] = 1;
        color="";
        #if (($2 in mapping_nodes) && ($3 in mapping_nodes)) {
        #    split(mapping_nodes[$2], x);
        #    split(mapping_nodes[$3], y);
        #    color=",color=green";
        #}
        if (($2"\t"$3 in mapping_nodes) && mapping_nodes[$2"\t"$3] != 0) {
            color=",color=green";
        }
        if ($4 == 0 && $7 == 0) {
            #printf("%s_%d_%s->%s_%d_%s[label=\"%d_%dR_%d\"%s];\n", $3, LEN_TBL[$3], POS_TBL[$3], $2, LEN_TBL[$2], POS_TBL[$2], $5 - $4 + 1, $10, mapping_nodes[$2"\t"$3], color);
            #printf("%s_%d_%s->%s_%d_%s[label=\"%d_%dR_%d\"%s];\n", $2, LEN_TBL[$2], POS_TBL[$2], $3, LEN_TBL[$3], POS_TBL[$3], $8 - $7 + 1, $10, mapping_nodes[$2"\t"$3], color);
        } else if ($4 == 0) {
            printf("%s_%d_%s->%s_%d_%s[label=\"%d_%d_%d\"%s];\n", $3, LEN_TBL[$3], POS_TBL[$3], $2, LEN_TBL[$2], POS_TBL[$2], $5 - $4 + 1, $10, mapping_nodes[$2"\t"$3], color);
        } else if ($7 == 0) {
            printf("%s_%d_%s->%s_%d_%s[label=\"%d_%d_%d\"%s];\n", $2, LEN_TBL[$2], POS_TBL[$2], $3, LEN_TBL[$3], POS_TBL[$3], $8 - $7 + 1, $10, mapping_nodes[$2"\t"$3], color);
        } else {
            printf("%s_%d->%s_%d[label=\"%d_%dF\"%s];\n", $3, LEN_TBL[$3], $2, LEN_TBL[$2], $5 - $4 + 1, $10, color);
            printf("%s_%d->%s_%d[label=\"%d_%dF\"%s];\n", $2, LEN_TBL[$2], $3, LEN_TBL[$3], $8 - $7 + 1, $10, color);
        }
    }
} END {
    for (i in LEN_TBL) {
        if (!(i in EDGE_NODE)) {
            printf("%s_%d;\n", i, LEN_TBL[i]);
        }
    }
    printf("}\n");
}
