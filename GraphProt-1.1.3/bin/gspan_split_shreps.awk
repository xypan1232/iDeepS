# print only graph parts of shreps (s); convert shreps to single graphs (t)
# call:
# awk -f convert_gspan.awk
BEGIN{state=start}
/^t/{state="t"}
/^u/{state="u"}
/^w/{state="w"}
/^s/{
	state="s";
	gsub(/^s/, "t");
	print
}
/^[^tuws]/&&state=="s"{print};

# gsub(/^V/, "v");
