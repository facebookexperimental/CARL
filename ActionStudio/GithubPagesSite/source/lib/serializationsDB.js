/**
 * IndexedDB persistence layer for CARL examples and definitions.
 * All serialized CARL objects are stored as binary blobs in a single
 * object store keyed by auto-incremented integer IDs.
 *
 * @module lib/serializationsDB
 */

/**
 * Wraps the browser's IndexedDB API to store and retrieve serialized CARL
 * objects (examples and definitions).
 */
export class SerializationsDB {
    static _DB_NAME = "carl_files";
    static _OBJECT_STORE_NAME = "serializations";
    _db;

    static async loadAsync() {
        const db = new SerializationsDB();
        db._db = await new Promise((resolve, reject) => {
            const request = window.indexedDB.open(SerializationsDB._DB_NAME, 3);
            request.onsuccess = evt => resolve(evt.target.result);
            request.onerror = () => {
                reject("Error fetching database: recordings will not be persisted in-browser");
            };
            request.onupgradeneeded = (evt) => {
                const db = evt.target.result;

                db.createObjectStore(SerializationsDB._OBJECT_STORE_NAME, { keyPath: "id", autoIncrement: true });
            };
        });

        db._db.onerror = (evt) => {
            console.error(`Database error: ${evt.target.error?.message}`);
        };

        return db;
    }

    static async resetAsync() {
        await new Promise((resolve, reject) => {
            const request = window.indexedDB.deleteDatabase(SerializationsDB._DB_NAME);
            request.onerror = () => reject("Error deleting database: we are in an unknown state");
            request.onsuccess = () => resolve();
        });
    }

    async addAsync(obj) {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readwrite")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        return await new Promise((resolve, reject) => {
            const request = os.add(obj);
            request.onerror = () => reject("Error storing CARL object");
            request.onsuccess = () => resolve(request.result);
        });
    }

    async fetchAllAsync() {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readonly")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        return await new Promise((resolve, reject) => {
            const request = os.getAll();
            request.onerror = () => reject("Error loading CARL objects");
            request.onsuccess = (evt) => {
                resolve(evt.target.result);
            };
        });
    }

    async updateAsync(id, updates) {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readwrite")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        let fetched = await new Promise((resolve, reject) => {
            const request = os.get(id);
            request.onerror = () => reject("Error in get() request");
            request.onsuccess = (evt) => {
                resolve(evt.target.result);
            };
        });
        fetched = { ...fetched, ...updates };

        return await new Promise((resolve, reject) => {
            const request = os.put(fetched);
            request.onerror = () => reject("Error in put() request");
            request.onsuccess = () => resolve();
        });
    }

    async deleteAsync(id) {
        const transaction = this._db.transaction(SerializationsDB._OBJECT_STORE_NAME, "readwrite")
        const os = transaction.objectStore(SerializationsDB._OBJECT_STORE_NAME);

        return await new Promise((resolve, reject) => {
            const request = os.delete(id);
            request.onerror = () => reject("Error deleting CARL objects");
            request.onsuccess = () => resolve();
        });
    }
}
