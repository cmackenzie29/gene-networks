function analyze(filename, genes, genenames)
	data = read("data/articles/$(filename)", String) |> 
				uppercase |> 
				x -> split(x, "\n\n\n") |>
				x -> map(y -> join(split(y, "\n\n")[2:end-1], " "), x) |>
				x -> map(y -> replace(y, " \n"=>" "), x) |>
				x -> map(y -> replace(y, "\n "=>" "), x) |>
				x -> map(y -> replace(y, "\n"=>" "), x) |>
				x -> filter(y -> length(y) >= 200, x)

	gene_mask = zeros(Bool, length(genes))
	genename_mask = zeros(Bool, length(genes))
	for g in 1:length(genes)
		for d in 1:length(data)
			if occursin(genes[g], data[d])
				gene_mask[g] = true
				break
			end
		end
	end
	for g in 1:length(genes)
		for d in 1:length(data)
			if occursin(genenames[g], data[d])
				genename_mask[g] = true
				break
			end
		end
	end

	hits = zeros(Int, length(genes))

	# Check instances with both
	for g in collect(1:length(genes))[gene_mask .& genename_mask]
		hits[g] += (map(y -> occursin(" GENE "*genes[g]*" ", y), data) .|| 
			map(y -> occursin(" GENE, "*genes[g]*" ", y), data) .|| 
			map(y -> occursin(" "*genes[g]*" GENE ", y), data) .||
			map(y -> occursin(" "*genes[g]*" EXPRESSION ", y), data) .||
			map(y -> occursin(" "*genenames[g]*" ", y), data)) |> sum
	end

	for g in collect(1:length(genes))[gene_mask .& .~genename_mask]
		hits[g] += (map(y -> occursin(" GENE "*genes[g]*" ", y), data) .|| 
			map(y -> occursin(" GENE, "*genes[g]*" ", y), data) .|| 
			map(y -> occursin(" "*genes[g]*" GENE ", y), data) .||
			map(y -> occursin(" "*genes[g]*" EXPRESSION ", y), data)) |> sum
	end

	for g in collect(1:length(genes))[.~gene_mask .& genename_mask]
		hits[g] += map(y -> occursin(" "*genenames[g]*" ", y), data) |> sum
	end

	return hits, length(data)
end

function analyze_set(disease_name)
	if map(x -> isfile("data/articles/pubmed_$(disease_name)_$(x).txt"), collect(2000:2025)) |> all

		genes = read("data/HuRI_genenames.tsv", String) |> 
					x -> split(x, r"[\n\t]") |> 
					unique |> 
					x -> filter(y -> y != "", x) |> 
					sort

		genenames = read("data/ensembl.tsv", String) |>
					x -> split(x, r"[\n\t]") |>
					x -> Dict(zip(x[4:3:end], x[5:3:end])) |>
					x -> [uppercase(x[g]) for g in genes]

		hits_total = zeros(Int, length(genes))
		n_articles_total = 0

		for year in 2000:2025
			println("Analyzing abstracts from $(year)...")
			hits, n_articles = analyze("pubmed_$(disease_name)_$(year).txt", genes, genenames)
			hits_total += hits
			n_articles_total += n_articles
		end

		output = ""
		for i in collect(1:length(genes))[hits .> 0]
			output = output * genes[i] * "\t" * "$(hits[i])" * "\n"
		end
		
		file = open("data/results/pubmed_$(disease_name)_results.txt", "w")
		write(file, output)
		close(file)

		println("Finished analysis of \"$(disease_name)\" with $(n_articles_total) articles")
	else
		println("Missing files with name data/articles/pubmed_$(disease_name)_YYYY.txt")
	end
end

print("\n\n\nEnter search term as it appears in the .txt files: ")
analyze_set(readline())


