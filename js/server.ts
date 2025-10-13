import { watch } from 'fs';
import { resolve, extname } from 'path';

const PORT = 3000;

// Define types for HMR clients
interface HMRClient {
  send: (data: string) => void;
}

const clients = new Set<HMRClient>();

// MIME type mapping
const mimeTypes: Record<string, string> = {
  '.html': 'text/html',
  '.js': 'text/javascript',
  '.css': 'text/css',
  '.json': 'application/json',
  '.png': 'image/png',
  '.jpg': 'image/jpeg',
  '.gif': 'image/gif',
  '.svg': 'image/svg+xml',
  '.wasm': 'application/wasm',
  '.ico': 'image/x-icon'
};

// Inject HMR client script into HTML files
function injectHMRScript(html: string): string {
  const hmrScript = `
  <script>
    // Hot Module Reloading client
    const eventSource = new EventSource('/__hmr');
    eventSource.onmessage = (event) => {
      if (event.data === 'reload') {
        console.log('[HMR] Reloading page...');
        location.reload();
      }
    };
    eventSource.onerror = () => {
      console.log('[HMR] Connection lost, retrying...');
    };
  </script>
  `;
  return html.replace('</body>', `${hmrScript}</body>`);
}

// Create the HTTP server
const server = Bun.serve({
  port: PORT,
  async fetch(req: Request): Promise<Response> {
    const url = new URL(req.url);

    // SSE endpoint for HMR
    if (url.pathname === '/__hmr') {
      const stream = new ReadableStream({
        start(controller: ReadableStreamDefaultController<Uint8Array>) {
          const encoder = new TextEncoder();

          // Add client to the set
          const client: HMRClient = {
            send: (data: string) => {
              controller.enqueue(encoder.encode(`data: ${data}\n\n`));
            }
          };
          clients.add(client);

          // Send initial connection message
          client.send('connected');

          // Clean up on close
          req.signal.addEventListener('abort', () => {
            clients.delete(client);
          });
        }
      });

      return new Response(stream, {
        headers: {
          'Content-Type': 'text/event-stream',
          'Cache-Control': 'no-cache',
          'Connection': 'keep-alive',
        },
      });
    }

    // Serve static files
    let pathname = url.pathname;
    if (pathname === '/') {
      pathname = '/index.html';
    }

    const filePath = resolve(import.meta.dir, '.' + pathname);
    const file = Bun.file(filePath);

    if (await file.exists()) {
      const ext = extname(pathname);
      const contentType = mimeTypes[ext] || 'application/octet-stream';

      // Inject HMR script into HTML files
      if (ext === '.html') {
        const html = await file.text();
        const modifiedHtml = injectHMRScript(html);
        return new Response(modifiedHtml, {
          headers: { 'Content-Type': contentType },
        });
      }

      return new Response(file, {
        headers: { 'Content-Type': contentType },
      });
    }

    return new Response('Not Found', { status: 404 });
  },
});

// Watch for file changes
const watchPaths: string[] = ['./src', './dist', './index.html'];

watchPaths.forEach((path: string) => {
  watch(path, { recursive: true }, (eventType: string, filename: string | null) => {
    if (filename) {
      console.log(`[HMR] File changed: ${filename}`);
      // Notify all connected clients
      clients.forEach((client: HMRClient) => {
        try {
          client.send('reload');
        } catch (err) {
          // Client disconnected, will be cleaned up
        }
      });
    }
  });
});

console.log(`Development server running at http://localhost:${PORT}`);
console.log(`Hot Module Reloading enabled`);
console.log(`Watching: ${watchPaths.join(', ')}`);
