/**
 * Shared utility functions used across the GithubPagesSite frontend.
 *
 * @module lib/utils
 */

/**
 * Triggers a browser download of a binary blob.
 *
 * @param {Uint8Array} jsBytes - The bytes to download.
 * @param {string} filename - The suggested filename for the download.
 */
export const downloadSerializedBytes = (jsBytes, filename) => {
    const buffer = jsBytes.buffer;
    const blob = new Blob([buffer], { type: "application/octet-stream" });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(link.href);
};

/**
 * Formats a duration in seconds as MM:SS.
 *
 * @param {number} seconds - Duration in seconds.
 * @returns {string} Formatted string, e.g. "01:23".
 */
export const formatDuration = (seconds) => {
    const mins = Math.floor(seconds / 60);
    const secs = Math.floor(seconds % 60);
    return `${mins.toString().padStart(2, '0')}:${secs.toString().padStart(2, '0')}`;
};
