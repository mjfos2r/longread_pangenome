from tqdm import tqdm
import time
import sys

def custom_write(text):
    # Use the tqdm.write method to ensure that the progress bar does not get disrupted
    tqdm.write(text)
    # Move the cursor up one line and clear the line
    sys.stdout.write('\033[F\033[K')

if __name__ == "__main__":
    # Example long-running task
    total_iterations = 100

    with tqdm(total=total_iterations) as pbar:
        for i in range(total_iterations):
            # Simulate a task
            time.sleep(0.1)
            # Update the progress bar
            pbar.update(1)
            # Update the line of text above the progress bar
            custom_write(f'Current iteration: {i+1}/{total_iterations}')

        # Final update to ensure the message stays after the loop finishes
        tqdm.write(f'Final iteration: {total_iterations}/{total_iterations}')