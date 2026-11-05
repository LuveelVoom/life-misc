package jls;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.prefs.Preferences;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;

// Sends the field being edited to Logic Life Search, which can use isotropic non-totalistic
// rules. The field goes through a temporary status file and LLS's own JLS status file import.
public class LlsRunDialog extends JDialog
{
    private static final long serialVersionUID = 1L;
    private static final String BACKTRACKER = "llsbt (backtracker)";

    private final ContentPane contentPane;
    private final Preferences prefs = Preferences.userNodeForPackage(LlsRunDialog.class);
    private final File llsDir = findLlsDir();

    private final JTextField ruleField = new JTextField(prefs.get("rule", ""), 24);
    private final JComboBox<String> solverBox = new JComboBox<String>();
    private final JSpinner solutionsSpinner = new JSpinner(new SpinnerNumberModel(1, 1, 1000000, 1));
    private final JSpinner parallelSpinner = new JSpinner(new SpinnerNumberModel(prefs.getInt("parallel", 4), 1, 64, 1));
    private final JCheckBox nonEmptyBox = new JCheckBox("Generation 0 must not be empty", prefs.getBoolean("nonEmpty", true));
    private final JCheckBox exactPeriodBox = new JCheckBox("Exclude lower periods (no shift or translation only)", prefs.getBoolean("exactPeriod", true));
    private final JTextField extraField = new JTextField(prefs.get("extra", ""), 24);
    private final JButton runButton = new JButton("Run");
    private final JButton stopButton = new JButton("Stop");
    private final JTextArea logArea = new JTextArea(24, 72);

    // the run in progress, used on the Swing thread except for stop()
    private final List<Process> processes = new ArrayList<Process>();
    private final List<File> outputCopies = new ArrayList<File>();
    private File statusFile;
    private File outputFile;
    private long startTime;
    private int copiesRunning;
    private int winner;
    private boolean running;
    private boolean finishing;

    public LlsRunDialog(JFrame owner, ContentPane contentPane)
    {
        super(owner, "Run in LLS", false);
        this.contentPane = contentPane;

        JPanel fields = new JPanel(new GridLayout(0, 2, 5, 3));
        fields.add(new JLabel("INT rule (blank: the rule from Properties)"));
        fields.add(ruleField);
        fields.add(new JLabel("SAT solver"));
        fields.add(solverBox);
        fields.add(new JLabel("Solutions to find"));
        fields.add(solutionsSpinner);
        fields.add(new JLabel("Parallel SAT copies, or backtracker threads"));
        fields.add(parallelSpinner);
        fields.add(new JLabel("Extra LLS arguments"));
        fields.add(extraField);
        fields.add(nonEmptyBox);
        fields.add(exactPeriodBox);

        JPanel buttons = new JPanel(new FlowLayout(FlowLayout.LEFT));
        buttons.add(runButton);
        buttons.add(stopButton);

        JPanel top = new JPanel(new BorderLayout());
        top.add(fields, BorderLayout.CENTER);
        top.add(buttons, BorderLayout.SOUTH);

        logArea.setEditable(false);
        logArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));

        add(top, BorderLayout.NORTH);
        add(new JScrollPane(logArea), BorderLayout.CENTER);

        for (String solver : findSolvers())
        {
            solverBox.addItem(solver);
        }
        if (new File(llsDir.getParentFile(), "native/llsbt").isFile())
        {
            solverBox.addItem(BACKTRACKER);
        }
        solverBox.setSelectedItem(prefs.get("solver", "kissat"));
        stopButton.setEnabled(false);
        runButton.addActionListener(e -> run());
        stopButton.addActionListener(e -> stop());
        getRootPane().setDefaultButton(runButton);

        // the solver must not outlive JLS
        Runtime.getRuntime().addShutdownHook(new Thread(this::stop));

        pack();
        setLocationRelativeTo(owner);
    }

    // LLS lives in the lls folder next to the jar
    private static File findLlsDir()
    {
        try
        {
            File codeSource = new File(LlsRunDialog.class.getProtectionDomain().getCodeSource().getLocation().toURI());
            return new File(codeSource.getParentFile(), "lls");
        }
        catch (URISyntaxException e)
        {
            return new File("lls");
        }
    }

    private List<String> findSolvers()
    {
        List<String> solvers = new ArrayList<String>();
        String[] names = new File(llsDir, "solvers").list();
        if (null != names)
        {
            Arrays.sort(names);
            for (String name : names)
            {
                if (!name.startsWith("README") && !name.startsWith("."))
                {
                    solvers.add(name);
                }
            }
        }
        return solvers;
    }

    // the LLS command line, apart from the status file name
    private List<String> makeArguments()
    {
        List<String> args = new ArrayList<String>();
        String rule = ruleField.getText().trim();
        if (rule.length() > 0)
        {
            args.add("--rule");
            args.add(rule);
        }
        args.add("-n");
        args.add(solutionsSpinner.getValue().toString());
        if (nonEmptyBox.isSelected())
        {
            args.add("--population");
            args.add(">0");
        }
        Properties properties = contentPane.getProperties();
        if (exactPeriodBox.isSelected() && properties.isTileGen()
                && (0 == properties.getTileGenShiftRight()) && (0 == properties.getTileGenShiftDown())
                && (Properties.XLAT_NONE == properties.getTranslation()))
        {
            int period = properties.getGenerations();
            for (int q = 2; q <= period; ++q)
            {
                if ((0 == period % q) && isPrime(q))
                {
                    args.add("--force_change");
                    args.add("0");
                    args.add(Integer.toString(period / q));
                }
            }
        }
        String extra = extraField.getText().trim();
        if (extra.length() > 0)
        {
            args.addAll(Arrays.asList(extra.split("\\s+")));
        }
        return args;
    }

    private static boolean isPrime(int n)
    {
        for (int d = 2; d * d <= n; ++d)
        {
            if (0 == n % d)
            {
                return false;
            }
        }
        return n > 1;
    }

    private void run()
    {
        if (running)
        {
            return;
        }
        if (!contentPane.isEditing())
        {
            JOptionPane.showMessageDialog(this, "LLS takes the field being edited. Reset the JLS search first.",
                    "Run in LLS", JOptionPane.INFORMATION_MESSAGE);
            return;
        }
        if (!new File(llsDir, "lls").isFile())
        {
            JOptionPane.showMessageDialog(this, "LLS not found at " + llsDir, "Run in LLS", JOptionPane.ERROR_MESSAGE);
            return;
        }
        String solver = (String) solverBox.getSelectedItem();
        int parallel = (Integer) parallelSpinner.getValue();
        prefs.put("rule", ruleField.getText());
        prefs.put("extra", extraField.getText());
        prefs.putBoolean("nonEmpty", nonEmptyBox.isSelected());
        prefs.putBoolean("exactPeriod", exactPeriodBox.isSelected());
        prefs.putInt("parallel", parallel);
        if (null != solver)
        {
            prefs.put("solver", solver);
        }

        processes.clear();
        outputCopies.clear();
        outputFile = null;
        try
        {
            statusFile = File.createTempFile("jls-lls-", ".jdf");
            contentPane.saveFieldForLls(statusFile.getPath());
            List<String> args = makeArguments();
            startTime = System.currentTimeMillis();
            if (BACKTRACKER.equals(solver))
            {
                logArea.setText("backtrack.py <field> --threads " + parallel + " " + String.join(" ", args) + "\n\n");
                List<String> command = new ArrayList<String>(Arrays.asList("python3", "-u", "backtrack.py", statusFile.getPath(),
                        "--threads", Integer.toString(parallel)));
                command.addAll(args);
                startCopy(command, llsDir.getParentFile());
            }
            else
            {
                logArea.setText("lls <field> --solver " + solver + " " + String.join(" ", args) + "\n");
                logArea.append(parallel > 1 ? parallel + " copies with different seeds, the first answer wins. Showing copy 1.\n\n" : "\n");
                int outputOption = Math.max(args.indexOf("-o"), args.indexOf("--output_file_name"));
                if ((parallel > 1) && (outputOption >= 0) && (outputOption + 1 < args.size()))
                {
                    outputFile = new File(args.get(outputOption + 1));
                }
                for (int copy = 0; copy < parallel; ++copy)
                {
                    List<String> command = new ArrayList<String>(Arrays.asList("python3", "-u", "lls", statusFile.getPath()));
                    if (null != solver)
                    {
                        command.add("--solver");
                        command.add(solver);
                    }
                    if (copy > 0)
                    {
                        // vary the seed, and let every other copy prefer dead cells
                        command.add("--parameters=--seed=" + copy + ((0 != copy % 2) ? " --phase=0" : ""));
                    }
                    List<String> copyArgs = new ArrayList<String>(args);
                    if (null != outputFile)
                    {
                        // copies must not mix their solutions in one file
                        File copyFile = new File(outputFile.getPath() + ".copy" + copy);
                        outputCopies.add(copyFile);
                        copyArgs.set(outputOption + 1, copyFile.getPath());
                    }
                    command.addAll(copyArgs);
                    startCopy(command, llsDir);
                }
            }
        }
        catch (IOException e)
        {
            logArea.append("Could not start LLS: " + e.getMessage() + "\n");
            stop();
            if (processes.isEmpty())
            {
                if (null != statusFile)
                {
                    statusFile.delete();
                }
                return;
            }
        }
        running = true;
        finishing = false;
        winner = 0;
        copiesRunning = processes.size();
        runButton.setEnabled(false);
        stopButton.setEnabled(true);
    }

    private void startCopy(List<String> command, File directory) throws IOException
    {
        final Process process = new ProcessBuilder(command).directory(directory).redirectErrorStream(true).start();
        final int copy;
        synchronized (processes)
        {
            copy = processes.size();
            processes.add(process);
        }
        Thread reader = new Thread(() -> readOutput(process, copy), "LLS output " + copy);
        reader.setDaemon(true);
        reader.start();
    }

    private void readOutput(Process process, int copy)
    {
        StringBuilder text = new StringBuilder();
        try (BufferedReader output = new BufferedReader(new InputStreamReader(process.getInputStream())))
        {
            String line;
            while (null != (line = output.readLine()))
            {
                text.append(line).append('\n');
                if (0 == copy)
                {
                    final String shown = line;
                    SwingUtilities.invokeLater(() -> logArea.append(shown + "\n"));
                }
            }
            process.waitFor();
        }
        catch (IOException | InterruptedException e)
        {
            // the output ends when the process is killed
        }
        final int exitValue = process.isAlive() ? -1 : process.exitValue();
        SwingUtilities.invokeLater(() -> copyFinished(copy, exitValue, text.toString()));
    }

    private void copyFinished(int copy, int exitValue, String text)
    {
        --copiesRunning;
        // LLS prints "Error" when the solver fails, and still exits with 0
        boolean answered = (0 == exitValue) && !text.contains("Error");
        if (answered && !finishing)
        {
            finishing = true;
            winner = copy;
            if (0 != copy)
            {
                logArea.append("\n[Copy " + (copy + 1) + " answered first]\n" + text);
            }
            stop();
        }
        if (0 == copiesRunning)
        {
            allCopiesFinished(exitValue);
        }
    }

    private void allCopiesFinished(int exitValue)
    {
        double seconds = (System.currentTimeMillis() - startTime) / 1000.0;
        statusFile.delete();
        // backtrack.py can't clean up after itself when it is killed
        for (String suffix : new String[] { ".cnf", ".watch", ".trans" })
        {
            new File(statusFile.getPath() + suffix).delete();
        }
        keepWinnersOutput();
        logArea.append(String.format("\n[LLS finished after %.1f s%s]\n", seconds,
                finishing ? "" : ", no answer, last exit code " + exitValue));
        logArea.setCaretPosition(logArea.getDocument().getLength());
        synchronized (processes)
        {
            processes.clear();
        }
        running = false;
        runButton.setEnabled(true);
        stopButton.setEnabled(false);
    }

    private void keepWinnersOutput()
    {
        try
        {
            if ((null != outputFile) && (winner < outputCopies.size()) && outputCopies.get(winner).isFile())
            {
                Files.write(outputFile.toPath(), Files.readAllBytes(outputCopies.get(winner).toPath()),
                        StandardOpenOption.CREATE, StandardOpenOption.APPEND);
            }
        }
        catch (IOException e)
        {
            logArea.append("Could not write " + outputFile + ": " + e.getMessage() + "\n");
        }
        for (File file : outputCopies)
        {
            file.delete();
        }
    }

    // kills the solver as well: killing Python alone would leave the solver running
    private void stop()
    {
        List<Process> running;
        synchronized (processes)
        {
            running = new ArrayList<Process>(processes);
        }
        for (Process process : running)
        {
            process.descendants().forEach(ProcessHandle::destroyForcibly);
            process.destroyForcibly();
        }
    }
}
