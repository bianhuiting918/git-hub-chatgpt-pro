$ErrorActionPreference = "Stop"

$ids = @("6EQE", "5XJH", "5YFE", "6ILW", "6EQD", "6EQF", "6EQG", "6EQH", "6QGC")
$outDir = "work/blind_work/01_system_setup/structures"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$rows = foreach ($id in $ids) {
  $url = "https://files.rcsb.org/download/$id.pdb"
  $path = Join-Path $outDir "$id.pdb"
  try {
    Invoke-WebRequest -Uri $url -OutFile $path
    $item = Get-Item -LiteralPath $path
    [PSCustomObject]@{
      pdb = $id
      source_url = $url
      local_path = $path
      bytes = $item.Length
      status = "downloaded"
    }
  } catch {
    [PSCustomObject]@{
      pdb = $id
      source_url = $url
      local_path = $path
      bytes = 0
      status = "failed: $($_.Exception.Message)"
    }
  }
}

$rows | ConvertTo-Csv -NoTypeInformation
