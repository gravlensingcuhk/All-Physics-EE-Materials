import typst, sys
try:
    typst.compile("tutor.typ", output="tutor.pdf", root=".")
    print("COMPILED OK")
except typst.TypstError as e:
    print("ERRORS:\n", e)
    sys.exit(1)
