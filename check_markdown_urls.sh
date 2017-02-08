#! /bin/bash
# Checks status of external urls in all markdown files in repository.
# Supply path as argument. Usage: "check_markdown_urls.sh path".

# Find all markdown files and use IFS to overcome whitespace issues
markdown_files=()
while IFS=  read -r -d $'\0'; do
    markdown_files+=("$REPLY")
done < <(find "$1" -name '*.md' -print0)

# Iterate over each markdown file and check urls
for file in "${markdown_files[@]}"; do
	# Grep all urls from markdown file supplied as argument (imperfect grep)
	urls_to_check=($(grep -o -E '.{0,0}http.{0,200}' "$file" | cut -d ')' -f 1))

	# Iterate over urls and obtain HTTP status code
	for url in "${urls_to_check[@]}"; do
		urlstatus=$(curl -o /dev/null --silent --head --write-out '%{http_code}' "$url" )
	    
	    # Output those urls that are not OK or redirected
	    if [ $urlstatus != '200' ] && [ $urlstatus != '303' ]; then
			echo "$urlstatus HTTP ERROR for [$url] in [$file]"
		fi
	done
done