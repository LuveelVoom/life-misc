tell application "Terminal"
  repeat with w in windows
    repeat with t in tabs of w
      if tty of t is "/dev/ttys001" then return history of t
    end repeat
  end repeat
end tell
