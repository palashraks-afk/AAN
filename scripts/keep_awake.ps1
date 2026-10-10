# Keeps Windows from going to sleep while the AAN pipeline runs, then stops by itself.
#
# It asks Windows (SetThreadExecutionState) not to idle-sleep for as long as this process is alive. It changes no power
# settings. It exits when the pipeline log says "pipeline finished", or after 30 hours. A closed lid still sleeps.
# Stop it any time by ending the powershell process, or by deleting nothing: it checks the log every minute.
#
#   powershell -NoProfile -WindowStyle Hidden -File scripts\keep_awake.ps1

Add-Type -Namespace Win32 -Name Power -MemberDefinition @'
[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);
'@

$ES_CONTINUOUS = 0x80000000
$ES_SYSTEM_REQUIRED = 0x00000001
$log = "D:\AAN_data\derived\pipeline.log"
$deadline = (Get-Date).AddHours(30)

while ((Get-Date) -lt $deadline) {
    [Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS -bor $ES_SYSTEM_REQUIRED) | Out-Null
    if ((Test-Path $log) -and ((Get-Content $log -Tail 3 -ErrorAction SilentlyContinue) -match "pipeline finished")) { break }
    Start-Sleep -Seconds 60
}
[Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS) | Out-Null
