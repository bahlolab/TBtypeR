import java.nio.file.Path

class Helpers {

    private static Class Nextflow = Class.forName('nextflow.Nextflow')

    static Path path(String filename) {
        Nextflow.file(filename, checkIfExists: true)
    }

    private static ArrayList<String> read_lines(Path path) {
        ArrayList<String> lines = path.toFile().readLines()
        if (lines.size() == 0) {
            Nextflow.error("Error: File '$path' is empty")
        }
        lines
    }

    /*
    required: list of required column names, fail if supplied and not present in file
    names: list of column names, will ignore header
    na_as_null: replace the string 'NA' with null (defualt true)
     */
    private static ArrayList<Map> read_delim(Map par = [:], Path path, String delim) {
        par.na_as_null = par.na_as_null ?: true

        ArrayList<String> lines = read_lines(path)
        ArrayList<String> names

        if (par.names) {
            names = par.names as ArrayList<String>
        } else {
            names = lines[0].split(delim)
            lines = lines.drop(1)
            if (lines.size() == 0) {
                Nextflow.error("Error: File '$path' contains header only")
            }
        }
        if (par.required) {
            if (!names.containsAll(par.required as ArrayList<String>)) {
                String missing = (par.required - names).join(', ')
                Nextflow.error("Error: Required filed(s) $missing missing from '$path'")
            }

        }
        lines.collect {
            ArrayList<String> fields = it.split(delim)
            if (fields.size() != names.size()) {
                Nextflow.error("Error: Line '$it' in '$path' has ${fields.size()} fields instead of ${names.size()}")
            }
            if (par.na_as_null) {
                fields = fields.collect { it == 'NA' ? null: it }
            }
            [names, fields].transpose().collectEntries { k, v -> [(k): v] }
        }
    }

    static ArrayList<Map> read_tsv(Map par = [:], Path path) {
        read_delim(par, path , '\t')
    }

    static ArrayList<Map> read_tsv(Map par = [:], String filename) {
        read_delim(par, path(filename) , '\t')
    }

    static ArrayList<Map> read_csv(Map par = [:], Path path) {
        read_delim(par, path , ',')
    }

    static ArrayList<Map> read_csv(Map par = [:], String filename) {
        read_delim(par, path(filename) , ',')
    }
}