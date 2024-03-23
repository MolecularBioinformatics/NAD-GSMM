# make sure you have inkscape installed, and make sure you're executing from the "images" folder

$directories = @(".\supplements", ".\parp", ".\")

# supplemental files
foreach ($directory in $directories) {
    $suppFiles = Get-ChildItem -Path $directory -Filter figure*.svg
    
    foreach ($file in $suppFiles) {
        $outputFilePath = Join-Path -Path $directory -ChildPath ($file.BaseName + ".pdf")
        & inkscape $file.FullName --export-filename=$outputFilePath
    }
}
