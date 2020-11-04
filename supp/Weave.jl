using Weave

filename = normpath("docs", "Table.jmd")

# Julia markdown to HTML
weave(filename; doctype = "md2html", out_path = :pwd)

# Julia markdown to PDF
weave(filename; doctype = "md2pdf", out_path = :pwd)

# Julia markdown to Pandoc markdown
weave(filename; doctype = "pandoc", out_path = :pwd)
