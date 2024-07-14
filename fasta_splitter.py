import tkinter as tk
from tkinter import filedialog, messagebox

def split_sequence(sequence, section_length=800):
    sections = []
    for i in range(0, len(sequence), section_length):
        section = sequence[i:i+section_length]
        sections.append(section)
    return sections

def label_sections(sequence, sections):
    labeled_sections = []
    for i, section in enumerate(sections):
        start_range = i * 800 + 1
        end_range = min((i + 1) * 800, len(sequence))
        labeled_section = f"Range: {start_range}-{end_range}\n{section}"
        labeled_sections.append(labeled_section)
    return labeled_sections

def browse_output_directory():
    output_directory = filedialog.askdirectory()
    output_directory_entry.delete(0, tk.END)
    output_directory_entry.insert(0, output_directory)

def split_and_save():
    sequence = sequence_entry.get().upper()
    output_directory = output_directory_entry.get()

    if not sequence:
        messagebox.showerror("Error", "Please enter a sequence")
        return

    sections = split_sequence(sequence)
    labeled_sections = label_sections(sequence, sections)

    output_filename = filedialog.asksaveasfilename(initialdir=output_directory, title="Save As", defaultextension=".txt", filetypes=[("Text files", "*.txt")])

    if output_filename:
        with open(output_filename, "w") as f:
            for labeled_section in labeled_sections:
                f.write(labeled_section + "\n\n")

        show_completion_message(output_filename)


def show_completion_message(output_directory):
    messagebox.showinfo("Success", f"Sequence split and saved to {output_directory}")

# Create main window
root = tk.Tk()
root.title("Amino Acid Sequence Splitter")

# Create input widgets
sequence_label = tk.Label(root, text="Enter the amino acid sequence:")
sequence_label.grid(row=0, column=0, sticky="w")
sequence_entry = tk.Entry(root, width=50)
sequence_entry.grid(row=0, column=1, padx=5, pady=5)

output_directory_label = tk.Label(root, text="Select output directory:")
output_directory_label.grid(row=1, column=0, sticky="w")
output_directory_entry = tk.Entry(root, width=40)
output_directory_entry.grid(row=1, column=1, padx=5, pady=5)
browse_button = tk.Button(root, text="Browse", command=browse_output_directory)
browse_button.grid(row=1, column=2, padx=5)

# Create button to split and save
split_button = tk.Button(root, text="Split and Save", command=split_and_save)
split_button.grid(row=2, column=0, columnspan=2, pady=10)

# Start GUI main loop
root.mainloop()
