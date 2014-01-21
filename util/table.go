package util

import (
	"errors"
	"fmt"
	"html/template"
	"io"
	"log"
	"os"
)

type Table struct {
	Title                  string
	ColHeaders, RowHeaders []string
	Data                   map[string][][]float64
}

func WriteTablesFile(tables []Table, filePath string) (err error) {
	file, err := os.Create(filePath)
	if err != nil {
		log.Println("opening file:", err)
		return
	}
	defer file.Close()

	return writeTablesHTML(tables, file)
}

func tableSanityCheck(table *Table) error {
	if table == nil {
		return errors.New("nil data table")
	}

	cols := len(table.ColHeaders)
	rows := len(table.RowHeaders)

	for _, dataSet := range table.Data {
		if actualRows := len(dataSet); actualRows != rows {
			return errors.New(fmt.Sprintf("inconsistent row counts: %v headers, %v rows", rows, actualRows))
		}
		for _, row := range dataSet {
			if len(row) != cols {
				return errors.New("inconsistent col counts")
			}
		}
	}

	return nil
}

func writeTablesHTML(tables []Table, output io.Writer) (err error) {
	for t := range tables {
		err = tableSanityCheck(&tables[t])
		if err != nil {
			return
		}
	}

	// Define a template.
	const document = `
<!DOCTYPE html>
<html>
<head>
    <style type="text/css">
        .results
        {
            font-family:"Trebuchet MS", Arial, Helvetica, sans-serif;
            width:100%;
            border-collapse:collapse;
        }
        .results td, .results th
        {
            font-size:1em;
            border:1px solid #98bf21;
            padding:3px 7px 2px 7px;
        }
        .results th
        {
            font-size:1.1em;
            text-align:left;
            padding-top:5px;
            padding-bottom:4px;
            background-color:#A7C942;
            color:#ffffff;
        }
        .results tr.alt td
        {
            color:#000000;
            background-color:#EAF2D3;
        }
        caption {
            text-align: left;
        }
    </style>
</head>
<body>
{{range $table := .}}
	<h2>{{.Title}}</h2>
	{{range $dataTitle, $data := $table.Data}}
	<table class="results">
	  <caption>{{$table.Title}} - {{$dataTitle}}</caption>
	  <tr>
	  	<th></th>
		{{range $table.ColHeaders}}<th>{{.}}</th>{{end}}
	  </tr>
	  {{range $index, $element := $data}}
	  <tr {{if eq (mod $index 2) 1}}class="alt"{{end}}>
		<th>{{index $table.RowHeaders $index}}</th>
		{{range $element}}<td>{{.}}</td>{{end}}
	  </tr>
	  {{end}}
	</table>
	{{end}}
{{end}}
</body>
</html>
`
	funcMap := template.FuncMap{
		"mod": func(a, b int) int { return a % b },
	}
	tDocument := template.Must(template.New("document").Funcs(funcMap).Parse(document))

	err = tDocument.Execute(output, tables)
	if err != nil {
		log.Println("executing template:", err)
	}

	return
}
