package util

import (
	"bytes"
	"os"
	"testing"
)

var tables = []Table{
	Table{
		Title:      "Test",
		RowHeaders: []string{"row1", "row2", "row3"},
		ColHeaders: []string{"col1", "col2"},
		Data: map[string][][]float64{
			"TestDataSet1": [][]float64{
				{0.0, 1.0},
				{2.0, 3.0},
				{4.0, 5.0},
			},
			"TestDataSet2": [][]float64{
				{1, 2},
				{1, 2},
				{1, 2},
			},
		},
	},
}

func TestTable(t *testing.T) {
	var b bytes.Buffer
	err := writeTablesHTML(tables, &b)
	if err != nil {
		t.Fatalf("error writing html %q", err)
	}

	t.Log(b.String())
}

func TestTablesFile(t *testing.T) {
	var fileName = "testTables.html"
	defer os.Remove("testTables.html")

	WriteTablesFile(tables, fileName)

	if file, err := os.OpenFile(fileName, os.O_RDONLY, os.ModePerm); err == nil {
		defer file.Close()

		if stat, err := file.Stat(); err == nil {
			if stat.Size() == 0 {
				t.Fatal("no text written")
			} else {
				t.Logf("Generated File is %v Bytes", stat.Size())
			}
		} else {
			t.Fatal(err)
		}
	} else {
		t.Fatal(err)
	}
}
